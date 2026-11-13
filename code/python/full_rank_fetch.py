#!/usr/bin/env python3
"""W1 full-ranking data: for every cohort molecule fetch from openFDA drug/event
  - the search string chosen by the pipeline's own field ladder (faers_search) and its report total
  - the counts of its top 1,000 reactions (patient.reaction.reactionmeddrapt.exact)
and the global background (total reports, top 1,000 reactions).
Writes data/raw/openfda_full/profiles.jsonl (one line per molecule, appended, resumable) and background.json.
One process, one call site. Read-only: touches no database."""
import csv, json, os, sys, time, threading, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", "")  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
sys.path.insert(0, os.path.join(ROOT, "openfda")); sys.path.insert(0, os.path.join(ROOT, "scripts"))
os.environ["FDA_API_KEY"] = subprocess.run([sys.executable, os.path.join(HERE, "fda_key.py")], capture_output=True, text=True).stdout.strip()
from openfda_client import OpenFDA                 # noqa: E402
from build_faers_signals import faers_search       # noqa: E402
OUT = os.path.join(PROJ, "data", "raw", "openfda_full")
RF = "patient.reaction.reactionmeddrapt.exact"
fda = OpenFDA(retries=8)
print("key:", "set" if fda.api_key else "none", flush=True)

def q(**kw):
    try:
        return fda.query("drug/event", **kw)
    except Exception as e:
        if "404" in str(e): return {}
        raise

bgp = os.path.join(OUT, "background.json")
if not os.path.exists(bgp):
    r = q(limit=1); g = q(count=RF, limit=1000)
    json.dump(dict(fetched=time.strftime("%Y-%m-%dT%H:%M:%S"), last_updated=r["meta"].get("last_updated"), N=r["meta"]["results"]["total"],
                   counts={x["term"]: x["count"] for x in g["results"]}), open(bgp, "w"))
print("background:", json.load(open(bgp))["N"], flush=True)

mols = [(r["ingredient_rxcui"], r["name"]) for r in csv.DictReader(open(os.path.join(PROJ, "data", "derived", "D2_cohort_molecules.csv")))]
cache = os.path.join(OUT, "profiles.jsonl"); done = set()
if os.path.exists(cache):
    for line in open(cache): done.add(json.loads(line)["ingredient_rxcui"])
todo = [m for m in mols if m[0] not in done]
print(f"molecules {len(mols)}, cached {len(done)}, to fetch {len(todo)}", flush=True)
lock = threading.Lock(); n = [0]; t0 = time.time()
def job(m):
    rx, name = m
    s, total = faers_search(fda, name)
    g = q(search=s, count=RF, limit=1000)
    rec = dict(ingredient_rxcui=rx, name=name, search=s, total=total, fetched=time.strftime("%Y-%m-%dT%H:%M:%S"),
               counts={x["term"]: x["count"] for x in (g.get("results") or [])})
    with lock:
        with open(cache, "a") as f: f.write(json.dumps(rec) + "\n")
        n[0] += 1
        if n[0] % 100 == 0: print(f"  {n[0]}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
with ThreadPoolExecutor(6) as ex:
    list(ex.map(job, todo))
print("done", sum(1 for _ in open(cache)), "profiles", flush=True)
