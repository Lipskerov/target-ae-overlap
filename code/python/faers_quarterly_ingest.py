#!/usr/bin/env python3
"""R10: ingest the FDA FAERS/AERS quarterly ASCII files (2004Q1 onward) into compact Parquet.

For each quarter: download the zip from fis.fda.gov, read DEMO / DRUG / REAC / INDI, keep the columns the
analysis needs, write one Parquet per table, delete the zip. Resumable (a finished quarter is skipped).
The server throttles each connection to about 60 KB/s, so quarters download in parallel (FAERS_WORKERS, default 14)
and are parsed one at a time.
The store lives OUTSIDE the repo (FAERS_STORE, default ~/faers_quarterly) because it is several GB; the
manifest (URL, bytes, sha256, row counts per quarter) is written into the project at
data/raw/faers_quarterly/MANIFEST.tsv so the source is reproducible.

Usage: faers_quarterly_ingest.py [first_quarter] [last_quarter]     e.g. 2004q1 2026q3
"""
import csv, hashlib, os, re, shutil, sys, threading, time, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly"))
MAN_DIR = os.path.join(PROJ, "data", "raw", "faers_quarterly"); os.makedirs(MAN_DIR, exist_ok=True)
MAN = os.path.join(MAN_DIR, "MANIFEST.tsv")
os.makedirs(os.path.join(STORE, "parquet"), exist_ok=True); TMP = os.path.join(STORE, "tmp")

def quarters(a, b):
    y, q = int(a[:4]), int(a[5]); out = []
    while (y, q) <= (int(b[:4]), int(b[5])):
        out.append(f"{y}q{q}"); q += 1
        if q == 5: y, q = y + 1, 1
    return out

def url_of(qt):   # legacy AERS naming up to 2012Q3
    return f"https://fis.fda.gov/content/Exports/{'aers' if qt <= '2012q3' else 'faers'}_ascii_{qt}.zip"

def find(names, prefix):
    c = [n for n in names if re.search(rf"(^|/){prefix}\d\dq\d(_new)?\.txt$", n, re.I)]   # 2018Q1 ships DEMO18Q1_new.txt
    return c[0] if c else None

WANT = {"demo": {"pid": ("primaryid", "isr"), "caseid": ("caseid", "case"), "ver": ("caseversion", "foll_seq"), "fda_dt": ("fda_dt",),
                 "init_fda_dt": ("init_fda_dt",), "occp_cod": ("occp_cod",), "rept_cod": ("rept_cod",)},
        "drug": {"pid": ("primaryid", "isr"), "drug_seq": ("drug_seq",), "role_cod": ("role_cod",), "drugname": ("drugname",), "prod_ai": ("prod_ai",)},
        "reac": {"pid": ("primaryid", "isr"), "pt": ("pt",)},
        "indi": {"pid": ("primaryid", "isr"), "drug_seq": ("indi_drug_seq", "drug_seq"), "indi_pt": ("indi_pt",)}}
INTS = {"pid", "caseid", "ver", "fda_dt", "init_fda_dt", "drug_seq"}

def parse(z, member, table):
    """Read one '$'-delimited member; returns (DataFrame with the wanted columns, lines skipped as malformed)."""
    with z.open(member) as f:
        total_lines = sum(1 for _ in f) - 1
    with z.open(member) as f:
        df = pd.read_csv(f, sep="$", dtype=str, quoting=3, encoding="latin-1", on_bad_lines="skip", index_col=False,
                         engine="c", low_memory=False, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    out = pd.DataFrame(index=df.index)
    for name, alts in WANT[table].items():
        src = next((a for a in alts if a in df.columns), None)
        col = df[src].str.strip() if src else pd.Series([""] * len(df), index=df.index)
        out[name] = pd.to_numeric(col, errors="coerce").astype("Int64") if name in INTS else col.str.upper()
    out = out[out["pid"].notna()]
    return out, total_lines - len(df), "isr" in df.columns

done = {}
if os.path.exists(MAN):
    for r in csv.DictReader(open(MAN), delimiter="\t"): done[r["quarter"]] = r
mlock, plock = threading.Lock(), threading.Lock()

def work(qt):
    outs = {t: os.path.join(STORE, "parquet", f"{t}_{qt}.parquet") for t in WANT}
    if qt in done and all(os.path.exists(p) for p in outs.values()): return
    t0 = time.time(); zp = os.path.join(TMP, f"{qt}.zip"); url = url_of(qt)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(4):
        try:
            h = hashlib.sha256(); size = 0
            with urllib.request.urlopen(req, timeout=180) as r, open(zp, "wb") as f:
                expect = int(r.headers.get("Content-Length") or 0)
                while True:
                    b = r.read(1 << 18)
                    if not b: break
                    f.write(b); h.update(b); size += len(b)
            if expect and size != expect: raise IOError(f"short read {size}/{expect}")
            break
        except Exception as e:
            print(qt, "download retry", attempt, e, flush=True); time.sleep(15)
    else:
        print(qt, "DOWNLOAD FAILED", flush=True); return
    td = time.time() - t0
    with plock:                                   # parse one quarter at a time (memory)
        try:
            with zipfile.ZipFile(zp) as z:
                names = z.namelist(); rows, bad, legacy = {}, {}, False
                for t in WANT:
                    m = find(names, t.upper())
                    if not m: raise IOError(f"missing member {t}")
                    df, bad[t], legacy = parse(z, m, t)
                    if t == "demo": df["legacy"] = legacy; df["quarter"] = qt
                    df.to_parquet(outs[t] + ".part", compression="zstd", index=False); os.replace(outs[t] + ".part", outs[t]); rows[t] = len(df)
        except Exception as e:
            print(qt, "PARSE FAILED", repr(e)[:300], flush=True); return
        finally:
            if os.path.exists(zp): os.remove(zp)
    with mlock:
        done[qt] = dict(quarter=qt, url=url, zip_bytes=size, zip_sha256=h.hexdigest(), legacy=int(legacy),
                        demo_rows=rows["demo"], drug_rows=rows["drug"], reac_rows=rows["reac"], indi_rows=rows["indi"],
                        skipped_lines=sum(bad.values()), ingested=time.strftime("%Y-%m-%dT%H:%M:%S"))
        with open(MAN, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(done[qt].keys()), delimiter="\t", extrasaction="ignore", restval=""); w.writeheader()
            for k in sorted(done): w.writerow(done[k])
        print(qt, f"zip {size/1e6:.0f}MB dl {td:.0f}s demo {rows['demo']:,} drug {rows['drug']:,} reac {rows['reac']:,} indi {rows['indi']:,} skipped {sum(bad.values())} [{len(done)} done]", flush=True)

first, last = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("2004q1", "2026q3")
shutil.rmtree(TMP, ignore_errors=True); os.makedirs(TMP)
todo = quarters(first, last)[::-1]               # newest (largest) first
with ThreadPoolExecutor(int(os.environ.get("FAERS_WORKERS", "14"))) as ex: list(ex.map(work, todo))
print("finished; quarters in manifest:", len(done), "of", len(quarters("2004q1", "2026q3")))
