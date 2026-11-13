#!/usr/bin/env python3
"""Scout for the co-reporting threat (PV review item 2): are same-target molecules named in the same
FAERS reports, and does the overlap survive among pairs that are rarely co-reported?
Seeded sample of same-target and different-target pairs from D10; one openFDA call per pair for the
number of reports matching BOTH molecules' search strings (strings from the full-ranking fetch).
Writes data/raw/openfda_full/coreport_sample.jsonl (resumable) and results/tables/coreport_sample.json.
Read-only; one process."""
import csv, gzip, json, os, random, subprocess, sys, time, threading
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", "")  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
sys.path.insert(0, os.path.join(ROOT, "openfda"))
os.environ["FDA_API_KEY"] = subprocess.run([sys.executable, os.path.join(HERE, "fda_key.py")], capture_output=True, text=True).stdout.strip()
from openfda_client import OpenFDA  # noqa: E402
fda = OpenFDA(retries=8)
N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SEED = 20261112
prof = {}
for line in open(os.path.join(PROJ, "data/raw/openfda_full/profiles.jsonl")):
    d = json.loads(line); prof[d["ingredient_rxcui"]] = d
same, diff = [], []
with gzip.open(os.path.join(PROJ, "data/derived/D10_pairs.csv.gz"), "rt") as f:
    for r in csv.DictReader(f):
        if r["a"] in prof and r["b"] in prof and prof[r["a"]]["total"] >= 100 and prof[r["b"]]["total"] >= 100:
            (same if r["same_target"] == "t" else diff).append((r["a"], r["b"], int(r["n_shared"]), int(r["n_shared_clean"]), r["shares_ind"] == "t"))
rng = random.Random(SEED)
sample = [("same",) + p for p in rng.sample(same, N_PER)] + [("diff",) + p for p in rng.sample(diff, N_PER)]
print(f"eligible pairs (both >=100 reports): same {len(same)}, diff {len(diff)}; sampling {N_PER} each", flush=True)
cache = os.path.join(PROJ, "data/raw/openfda_full/coreport_sample.jsonl"); done = {}
if os.path.exists(cache):
    for line in open(cache):
        d = json.loads(line); done[(d["a"], d["b"])] = d
lock = threading.Lock()
def job(s):
    g, a, b, n, nc, si = s
    if (a, b) in done: return
    q = f'({prof[a]["search"]}) AND ({prof[b]["search"]})'
    try:
        r = fda.query("drug/event", search=q, limit=1); both = r.get("meta", {}).get("results", {}).get("total", 0) or 0
    except Exception as e:
        if "404" in str(e): both = 0
        else: raise
    rec = dict(group=g, a=a, b=b, name_a=prof[a]["name"], name_b=prof[b]["name"], n_shared=n, n_shared_clean=nc, shares_ind=si,
               total_a=prof[a]["total"], total_b=prof[b]["total"], both=both)
    with lock:
        done[(a, b)] = rec
        with open(cache, "a") as f: f.write(json.dumps(rec) + "\n")
with ThreadPoolExecutor(6) as ex: list(ex.map(job, sample))
rows = [done[(s[1], s[2])] for s in sample]
for r in rows: r["share"] = r["both"] / min(r["total_a"], r["total_b"])
def mean(x): return round(sum(x) / len(x), 3) if x else None
def med(x): x = sorted(x); return round(x[len(x) // 2], 4) if x else None
out = dict(seed=SEED, n_per_group=N_PER, definition="co-report share = reports matching both molecules / the smaller molecule's report total", groups={}, strata=[])
for g in ("same", "diff"):
    rs = [r for r in rows if r["group"] == g]
    out["groups"][g] = dict(pairs=len(rs), median_share=med([r["share"] for r in rs]), mean_share=mean([r["share"] for r in rs]),
                            pct_share_ge_5pct=round(100 * sum(r["share"] >= 0.05 for r in rs) / len(rs), 1),
                            pct_zero_coreports=round(100 * sum(r["both"] == 0 for r in rs) / len(rs), 1), mean_shared=mean([r["n_shared"] for r in rs]))
for lo, hi in ((0, 0.001), (0.001, 0.01), (0.01, 0.05), (0.05, 2)):
    st = dict(share_from=lo, share_to=hi)
    for g in ("same", "diff"):
        rs = [r for r in rows if r["group"] == g and lo <= r["share"] < hi]
        st[g] = dict(pairs=len(rs), mean_shared=mean([r["n_shared"] for r in rs]), mean_shared_clean=mean([r["n_shared_clean"] for r in rs]))
    out["strata"].append(st)
json.dump(out, open(os.path.join(PROJ, "results/tables/coreport_sample.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
