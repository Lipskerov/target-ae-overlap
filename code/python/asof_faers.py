#!/usr/bin/env python3
"""As-of-date hindsight test for the DDO target piece.

Question: do the FAERS signals of a target's FIRST-approved molecule, as they stood on the day
each LATER molecule at that target was approved, predict the later molecule's CURRENT top-20
signals -- better than the as-of signals of first-approved molecules at OTHER targets?

Inputs (asof/*.tsv, exported by asof/00_export.sql from the same temp tables as the article):
  pairs.tsv   tgt, f_rx, f_name, f_ap, s_rx, s_name, s_ap
  firsts.tsv  tgt, rx, name, ap           (one first-approved molecule per cohort target)
  sig.tsv     rx, reaction_pt, kind       (current top-20 signals, raw and clean)
  ind.tsv     rx, indication              (current FAERS top indications, for cleaning)
  artefact.tsv reaction_pt                (the article's artefact term list)
Method: openFDA drug/event counts restricted to receivedate < s_ap; PRR/chi2 exactly as
scripts/compute_faers_prr.py (PRR>=2, a>=3, Yates chi2>=4, BLACKLIST), top-20 by PRR.
Baseline: for each successor, K first-approved molecules of OTHER targets approved before s_ap,
profiled as-of the same date (seeded sample). Writes asof/profiles.jsonl incrementally (resumable)
and asof/RESULT.txt. Read-only against the DB (reads nothing from it; uses the TSV exports).
"""
import csv, json, os, random, re, sys, time, threading
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from datetime import date, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("TRIALYTIC_ROOT", "")  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
sys.path.insert(0, os.path.join(ROOT, "openfda"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from openfda_client import OpenFDA                        # noqa: E402
from build_faers_signals import faers_search              # noqa: E402  (same field ladder as the stored profiles)
from compute_faers_prr import BLACKLIST                   # noqa: E402

K_BASE = int(os.environ.get("ASOF_K", "3"))
SEED = 20261021
FAERS_START = date(2004, 1, 1)   # openFDA drug/event coverage starts 2004 (checked live, see RESULT.txt)
MIN_REPORTS = (1, 100, 1000)     # strata for the first molecule's as-of report count

def read(name):
    with open(os.path.join(HERE, name), newline="") as f:
        return [r for r in csv.reader(f, delimiter="\t") if r and r[0]]

pairs = [dict(tgt=r[0], f_rx=r[1], f_name=r[2], f_ap=date.fromisoformat(r[3]), s_rx=r[4], s_name=r[5], s_ap=date.fromisoformat(r[6])) for r in read("pairs.tsv")]
firsts = [dict(tgt=r[0], rx=r[1], name=r[2], ap=date.fromisoformat(r[3])) for r in read("firsts.tsv")]
name_of = {f["rx"]: f["name"] for f in firsts}
sig = defaultdict(set)
for rx, rpt, kind in read("sig.tsv"):
    sig[(rx, kind)].add(rpt.upper())
own_ind = defaultdict(set)
for rx, ind in read("ind.tsv"):
    own_ind[rx].add(ind.upper())
artefact = {r[0].upper() for r in read("artefact.tsv")}
ART_RE = re.compile(r"(PRODUCT|OFF LABEL|DRUG INEFFECTIVE|DRUG INTERACTION|DRUG DOSE|DRUG ABUSE|DRUG DEPENDENCE|CONTRAINDICATED|INAPPROPRIATE|MEDICATION ERROR|NO ADVERSE|INTENTIONAL|WRONG|DOSE OMISSION|THERAPY |THERAPY$|TREATMENT|^DEATH$|PROGRESSION|METASTAS|COVID|^ILLNESS$|CONDITION AGGRAVATED|HOSPITALISATION|UNDERDOSE|OVERDOSE|EXPOSURE|INJECTION SITE|INFUSION SITE|APPLICATION SITE|IMPLANT SITE)")

fda = OpenFDA(retries=8)
calls = 0
def q(**kw):
    global calls
    calls += 1
    try:
        return fda.query("drug/event", **kw)
    except Exception as e:   # 404 = no matches
        if "404" in str(e): return {}
        raise
def cnt(field, search, limit=1000):
    r = q(search=search, count=field, limit=limit)
    return {x["term"].upper(): x["count"] for x in (r.get("results") or [])}
def total(search):
    r = q(search=search, limit=1)
    return (r.get("meta", {}).get("results", {}).get("total", 0) or 0)

# --- caches on disk ---------------------------------------------------------
CACHE = os.path.join(HERE, "profiles.jsonl")
prof = {}
if os.path.exists(CACHE):
    for line in open(CACHE):
        d = json.loads(line); prof[d["key"]] = d
LOCK = threading.Lock()
def save(d):
    with LOCK:
        prof[d["key"]] = d
        with open(CACHE, "a") as f: f.write(json.dumps(d) + "\n")

search_of = {}
SLOCKS = defaultdict(threading.Lock)
def search_for(rx):
  with SLOCKS[rx]:
    if rx in search_of: return search_of[rx]
    key = f"search|{rx}"
    if key in prof:
        search_of[rx] = (prof[key]["search"], prof[key]["total_now"]); return search_of[rx]
    s, t = faers_search(fda, name_of[rx]); global calls; calls += 4
    save(dict(key=key, rx=rx, name=name_of[rx], search=s, total_now=t)); search_of[rx] = (s, t); return s, t

def datestr(d): return d.strftime("%Y%m%d")
def before(d):  # receivedate strictly before d
    return f"receivedate:[19000101 TO {datestr(d - timedelta(days=1))}]"

BGLOCKS = defaultdict(threading.Lock)
def background(d):
    key = f"bg|{d.isoformat()}"
    with BGLOCKS[key]:
        if key in prof: return prof[key]
        N = total(before(d)); bg = cnt("patient.reaction.reactionmeddrapt.exact", before(d))
        rec = dict(key=key, N=N, bg=bg); save(rec); return rec

def profile(rx, d):
    """Top-20 PRR signals of molecule rx using only reports received before d."""
    key = f"p|{rx}|{d.isoformat()}"
    if key in prof: return prof[key]
    s, _ = search_for(rx)
    S = f"({s}) AND {before(d)}"
    ab = total(S)
    rec = dict(key=key, rx=rx, date=d.isoformat(), total=ab, top=[], top_clean=[], missing_bg=0)
    if ab >= 3:
        counts = cnt("patient.reaction.reactionmeddrapt.exact", S)
        bgrec = background(d); N, bg = bgrec["N"], bgrec["bg"]
        # Reactions outside the as-of global top-1000 have no background and are not scored --
        # the same limitation as the stored profiles (compute_faers_prr.py counts them as missing_bg).
        # Fetching them one by one cost ~40 calls per profile (first attempt, 2026-10-21) and was dropped.
        scored = []
        for r, a in counts.items():
            ac = bg.get(r)
            if not ac: rec["missing_bg"] += 1; continue
            b, c = ab - a, ac - a; dd = N - ab - c
            if b <= 0 or c <= 0 or dd <= 0: continue
            prr = (a / ab) / (c / (c + dd))
            exp = [ab * ac / N, ab * (N - ac) / N, (N - ab) * ac / N, (N - ab) * (N - ac) / N]
            chi2 = sum((abs(o - e) - 0.5) ** 2 / e for o, e in zip([a, b, c, dd], exp) if e > 0)
            if round(prr, 2) >= 2 and a >= 3 and round(chi2, 1) >= 4 and r not in BLACKLIST:
                scored.append((prr, a, r))
        scored.sort(key=lambda x: (-x[0], -x[1], x[2]))
        rec["n_signals"] = len(scored)
        rec["top"] = [r for _, _, r in scored[:20]]
        clean = [r for _, _, r in scored if not ART_RE.search(r) and r not in artefact and r not in own_ind[rx]]
        rec["top_clean"] = clean[:20]
    save(rec); return rec

# --- run ----------------------------------------------------------------------
rng = random.Random(SEED)
scored_pairs = [p for p in pairs if p["s_ap"] >= FAERS_START]
print(f"pairs {len(pairs)}, scoring {len(scored_pairs)} with successor approved >= {FAERS_START}", flush=True)
t0 = time.time(); rows = []
# work list first (deterministic sampling), then parallel fetch, then scoring from the cache
work, base_of = [], {}
for p in scored_pairs:
    cands = [f for f in firsts if f["tgt"] != p["tgt"] and f["rx"] != p["f_rx"] and f["ap"] < p["s_ap"]]
    base_of[(p["f_rx"], p["s_rx"])] = rng.sample(cands, min(K_BASE, len(cands)))
    work.append((p["f_rx"], p["s_ap"])); work += [(f["rx"], p["s_ap"]) for f in base_of[(p["f_rx"], p["s_rx"])]]
work = sorted(set(work), key=lambda w: (w[1], w[0]))
print(f"profiles to fetch: {len(work)} ({sum(1 for w in work if f'p|{w[0]}|{w[1].isoformat()}' in prof)} cached)", flush=True)
with ThreadPoolExecutor(int(os.environ.get("ASOF_THREADS", "6"))) as ex:   # search strings (4 probes each)
    list(ex.map(search_for, sorted({w[0] for w in work})))
print(f"search strings ready: {len(search_of)}, {calls} calls, {time.time()-t0:.0f}s", flush=True)
done = [0]
def job(w):
    profile(*w); done[0] += 1
    if done[0] % 50 == 0: print(f"  {done[0]}/{len(work)} profiles, {calls} calls, {time.time()-t0:.0f}s", flush=True)
with ThreadPoolExecutor(int(os.environ.get("ASOF_THREADS", "6"))) as ex:
    list(ex.map(job, work))
for i, p in enumerate(scored_pairs):
    own = profile(p["f_rx"], p["s_ap"])
    base = base_of[(p["f_rx"], p["s_rx"])]
    brec = [profile(f["rx"], p["s_ap"]) for f in base]
    def shared(rec, kind):
        top = set(rec["top"] if kind == "raw" else rec["top_clean"])
        return len(top & sig[(p["s_rx"], kind)])
    rows.append(dict(tgt=p["tgt"], f=p["f_name"], s=p["s_name"], s_ap=p["s_ap"].isoformat(), f_total=own["total"],
                     own_raw=shared(own, "raw"), own_clean=shared(own, "clean"),
                     base=[dict(rx=b["rx"], total=b["total"], raw=shared(b, "raw"), clean=shared(b, "clean")) for b in brec],
                     s_has_raw=len(sig[(p["s_rx"], "raw")]), s_has_clean=len(sig[(p["s_rx"], "clean")])))
json.dump(rows, open(os.path.join(HERE, "pairs_scored.json"), "w"), indent=0)

# --- summarise ----------------------------------------------------------------
def mean(xs): return sum(xs) / len(xs) if xs else float("nan")
out = []
out.append(f"as-of hindsight test, run {date.today().isoformat()}; openFDA calls {calls}; K_BASE={K_BASE}; seed {SEED}")
out.append(f"pairs total {len(pairs)}; successor approved before {FAERS_START}: {len(pairs)-len(scored_pairs)} (no FAERS as-of data, excluded); scored {len(scored_pairs)}")
out.append(f"background check: FAERS reports received before 2004-01-01 = {total(before(FAERS_START)):,}; before 2005-01-01 = {total(before(date(2005,1,1))):,}")
for thr in MIN_REPORTS:
    sel = [r for r in rows if r["f_total"] >= thr]
    bsel = [b for r in sel for b in r["base"] if b["total"] >= thr]
    out.append(f"\n[first molecule had >= {thr} FAERS reports as-of the successor's approval]  pairs={len(sel)} baseline-pairs={len(bsel)}")
    for kind in ("raw", "clean"):
        own = mean([r[f"own_{kind}"] for r in sel]); bas = mean([b[kind] for b in bsel])
        out.append(f"  {kind:5s}: same target as-of {own:.2f} | other-target firsts as-of {bas:.2f} | ratio {own/bas if bas else float('nan'):.1f}x"
                   f" | same-target share>=1: {100*mean([(r[f'own_{kind}']>=1) for r in sel]):.1f}% vs {100*mean([(b[kind]>=1) for b in bsel]):.1f}%")
# per-target view (>=100 reports stratum, raw)
sel = [r for r in rows if r["f_total"] >= 100]
per = defaultdict(list)
for r in sel: per[r["tgt"]].append(r["own_raw"])
bas_all = mean([b["raw"] for r in sel for b in r["base"] if b["total"] >= 100])
above = sum(mean(v) > bas_all for v in per.values())
out.append(f"\nper target (>=100 stratum, raw): {len(per)} targets, {above} above the baseline mean {bas_all:.2f}; median target mean {sorted(mean(v) for v in per.values())[len(per)//2]:.2f}")
out.append("\nexamples (raw, >=100 stratum), first -> successor (approved): shared as-of | baseline mean")
for r in sorted(sel, key=lambda r: -r["own_raw"])[:12]:
    out.append(f"  {r['tgt']:8s} {r['f']} -> {r['s']} ({r['s_ap']}): {r['own_raw']} | {mean([b['raw'] for b in r['base']]):.1f}  (first had {r['f_total']:,} reports then)")
open(os.path.join(HERE, "RESULT.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
