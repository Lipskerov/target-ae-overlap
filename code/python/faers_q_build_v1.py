#!/usr/bin/env python3
"""R10 step 2: from the quarterly Parquet store build the analysis database (SAP.md §1).

  cases     one row per case after deduplication (pid kept, caseid, initial receipt date, reporter type)
  drug_rx   (pid, rx, role_rank, fdc_only): each deduplicated report x RxNorm ingredient named in it
            role_rank 1 = primary suspect, 2 = secondary suspect, 3 = concomitant, 4 = interacting
            fdc_only  = the ingredient appears in that report only inside a fixed-dose-combination product
  reac      (pid, pt) distinct, deduplicated reports only
  indi      (pid, rx, indi_pt) indication attached to a mapped drug row

Writes ~/faers_quarterly/faers.duckdb (outside the repo) and, inside the project,
data/derived/faers_q/report_flow.json and name_map.csv.gz (every drug string -> ingredient(s), with row counts).
Read-only against everything else."""
import csv, gzip, json, os, re, sys, time
import duckdb

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", "")  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from build_drug_fda_crosswalk import normalize   # noqa: E402  the repository's canonical normalizer

STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly"))
PQ = os.path.join(STORE, "parquet"); OUT = os.path.join(PROJ, "data", "derived", "faers_q"); os.makedirs(OUT, exist_ok=True)
DB = os.path.join(STORE, "faers.duckdb")
if os.path.exists(DB): os.remove(DB)
con = duckdb.connect(DB)
con.execute(f"SET threads=6; SET memory_limit='5GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
flow = {}; t0 = time.time()
def one(q): return con.execute(q).fetchone()[0]
def step(msg): print(f"[{time.time()-t0:5.0f}s] {msg}", flush=True)

n_q = len([f for f in os.listdir(PQ) if f.startswith("demo_")]); flow["quarters"] = n_q
con.execute(f"CREATE VIEW demo_raw AS SELECT * FROM read_parquet('{PQ}/demo_*.parquet', union_by_name=true)")
flow["demo_rows_all_versions"] = one("SELECT count(*) FROM demo_raw")
flow["distinct_caseid"] = one("SELECT count(DISTINCT caseid) FROM demo_raw WHERE caseid IS NOT NULL")
flow["demo_rows_without_caseid"] = one("SELECT count(*) FROM demo_raw WHERE caseid IS NULL")
# Deduplicate: FAERS-era version beats legacy; then latest fda_dt; then highest pid. Initial receipt = earliest date seen.
con.execute("""CREATE TABLE cases AS
SELECT pid, caseid, occp_cod, legacy, quarter, fda_dt, init_dt FROM (
  SELECT *, min(coalesce(init_fda_dt, fda_dt)) OVER (PARTITION BY caseid) AS init_dt,
         row_number() OVER (PARTITION BY caseid ORDER BY legacy ASC, fda_dt DESC NULLS LAST, pid DESC) AS rn
  FROM demo_raw WHERE caseid IS NOT NULL) WHERE rn = 1""")
flow["cases_after_dedup"] = one("SELECT count(*) FROM cases")
flow["cases_legacy_era"] = one("SELECT count(*) FROM cases WHERE legacy")
flow["cases_without_receipt_date"] = one("SELECT count(*) FROM cases WHERE init_dt IS NULL")
step(f"cases {flow['cases_after_dedup']:,} of {flow['demo_rows_all_versions']:,} rows")

con.execute(f"""CREATE TABLE drug_d AS SELECT d.pid, d.drug_seq, d.role_cod, d.drugname, d.prod_ai
FROM read_parquet('{PQ}/drug_*.parquet', union_by_name=true) d SEMI JOIN cases c ON c.pid = d.pid""")
flow["drug_rows_in_dedup_cases"] = one("SELECT count(*) FROM drug_d")
step(f"drug rows {flow['drug_rows_in_dedup_cases']:,}")

# ---- name map: distinct strings -> ingredient(s) --------------------------------------------------
name2rx = {}
for r in csv.DictReader(open(os.path.join(PROJ, "data/raw/trialytic_db/R2_name_to_ingredient_map.csv"))):
    if r["ingredient_rxcui"]: name2rx[r["normalized_name"]] = r["ingredient_rxcui"]
SPLIT = re.compile(r"\\|/|;|\+|,|\bAND\b|\bWITH\b|\bW/", re.I)
cache = {}
def lookup(s):
    if s not in cache:
        cache[s] = name2rx.get(normalize(s)) if s else None
    return cache[s]
def resolve(prod_ai, drugname):
    """-> (tuple of rxcui, source). prod_ai is FDA's own active-ingredient field; '\\' separates components."""
    if prod_ai:
        parts = [p.strip() for p in prod_ai.split("\\") if p.strip()]
        hit = [lookup(p) for p in parts]
        if any(hit): return tuple(sorted({h for h in hit if h})), ("prod_ai_fdc" if len(parts) > 1 else "prod_ai")
    if drugname:
        h = lookup(drugname)
        if h: return (h,), "drugname"
        parts = [p.strip() for p in SPLIT.split(drugname) if p.strip()]
        if len(parts) > 1:
            hit = {lookup(p) for p in parts} - {None}
            if hit: return tuple(sorted(hit)), ("drugname_fdc" if len(parts) > 1 else "drugname")
    return (), "unmapped"

strings = con.execute("SELECT prod_ai, drugname, count(*) n FROM drug_d GROUP BY 1, 2").fetchall()
step(f"distinct (prod_ai, drugname) strings {len(strings):,}")
rows = []; src_rows = {}
with gzip.open(os.path.join(OUT, "name_map.csv.gz"), "wt", newline="") as f:
    w = csv.writer(f); w.writerow(["prod_ai", "drugname", "drug_rows", "ingredient_rxcuis", "source", "multi_component"])
    for pa, dn, n in strings:
        rx, src = resolve(pa or "", dn or "")
        multi = src.endswith("_fdc")
        src_rows[src] = src_rows.get(src, 0) + n
        if rx or n >= 50: w.writerow([pa, dn, n, "|".join(rx), src, int(multi)])   # rare unmapped strings are counted, not listed
        for x in rx: rows.append((pa, dn, x, multi))
flow["drug_rows_by_mapping_source"] = src_rows
con.execute("CREATE TABLE smap (prod_ai VARCHAR, drugname VARCHAR, rx VARCHAR, multi BOOLEAN)")
import pandas as pd
smap_df = pd.DataFrame(rows, columns=["prod_ai", "drugname", "rx", "multi"]); con.execute("INSERT INTO smap SELECT * FROM smap_df")
step(f"string map rows {len(rows):,}")

con.execute("""CREATE TABLE drug_rows_rx AS
SELECT d.pid, d.drug_seq, m.rx, m.multi,
       CASE d.role_cod WHEN 'PS' THEN 1 WHEN 'SS' THEN 2 WHEN 'C' THEN 3 WHEN 'I' THEN 4 ELSE 3 END AS role_rank
FROM drug_d d JOIN smap m ON m.prod_ai IS NOT DISTINCT FROM d.prod_ai AND m.drugname IS NOT DISTINCT FROM d.drugname""")
con.execute("""CREATE TABLE drug_rx AS
SELECT pid, rx, min(role_rank) AS role_rank, bool_and(multi) AS fdc_only FROM drug_rows_rx GROUP BY 1, 2""")
flow["report_x_ingredient_rows"] = one("SELECT count(*) FROM drug_rx")
flow["reports_with_a_mapped_drug"] = one("SELECT count(DISTINCT pid) FROM drug_rx")
flow["ingredients_seen"] = one("SELECT count(DISTINCT rx) FROM drug_rx")
step(f"drug_rx {flow['report_x_ingredient_rows']:,}")

con.execute(f"""CREATE TABLE reac AS SELECT DISTINCT r.pid, r.pt
FROM read_parquet('{PQ}/reac_*.parquet', union_by_name=true) r SEMI JOIN cases c ON c.pid = r.pid WHERE r.pt <> ''""")
flow["report_x_reaction_rows"] = one("SELECT count(*) FROM reac"); flow["distinct_pt"] = one("SELECT count(DISTINCT pt) FROM reac")
flow["reports_with_a_reaction"] = one("SELECT count(DISTINCT pid) FROM reac")
step(f"reac {flow['report_x_reaction_rows']:,}")

con.execute(f"""CREATE TABLE indi AS SELECT DISTINCT i.pid, x.rx, i.indi_pt
FROM read_parquet('{PQ}/indi_*.parquet', union_by_name=true) i
JOIN drug_rows_rx x ON x.pid = i.pid AND x.drug_seq = i.drug_seq WHERE i.indi_pt <> ''""")
flow["report_x_ingredient_x_indication_rows"] = one("SELECT count(*) FROM indi")
con.execute("DROP TABLE drug_d; DROP TABLE drug_rows_rx")
con.execute("CHECKPOINT"); con.close()
flow["built"] = time.strftime("%Y-%m-%dT%H:%M:%S"); flow["db_bytes"] = os.path.getsize(DB)
json.dump(flow, open(os.path.join(OUT, "report_flow.json"), "w"), indent=1)
print(json.dumps(flow, indent=1))
