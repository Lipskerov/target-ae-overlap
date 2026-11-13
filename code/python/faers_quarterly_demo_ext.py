#!/usr/bin/env python3
"""Second pass over the FDA quarterly files: DEMO fields needed for field-based duplicate detection and for the
case-characteristics table (age, sex, event date, reporter/occurrence country). Each zip is downloaded again, its
sha256 compared with data/raw/faers_quarterly/MANIFEST.tsv (same file as the first pass), and only the DEMO
member is parsed -> ~/faers_quarterly/parquet/demoext_<quarter>.parquet. Resumable."""
import csv, hashlib, os, re, shutil, sys, threading, time, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly")); TMP = os.path.join(STORE, "tmp2")
MAN = {r["quarter"]: r for r in csv.DictReader(open(os.path.join(PROJ, "data/raw/faers_quarterly/MANIFEST.tsv")), delimiter="\t")}
LOG = os.path.join(PROJ, "data/raw/faers_quarterly/MANIFEST_demoext.tsv")
WANT = {"pid": ("primaryid", "isr"), "caseid": ("caseid", "case"), "event_dt": ("event_dt",), "age": ("age",), "age_cod": ("age_cod",),
        "sex": ("sex", "gndr_cod"), "reporter_country": ("reporter_country",), "occr_country": ("occr_country",), "wt": ("wt",), "mfr_num": ("mfr_num",), "mfr_sndr": ("mfr_sndr",)}
done = {}
if os.path.exists(LOG):
    for r in csv.DictReader(open(LOG), delimiter="\t"): done[r["quarter"]] = r
lock, plock = threading.Lock(), threading.Lock()
def work(qt):
    out = os.path.join(STORE, "parquet", f"demoext_{qt}.parquet")
    if qt in done and os.path.exists(out): return
    zp = os.path.join(TMP, f"{qt}.zip"); req = urllib.request.Request(MAN[qt]["url"], headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(4):
        try:
            h = hashlib.sha256(); size = 0
            with urllib.request.urlopen(req, timeout=180) as r, open(zp, "wb") as f:
                expect = int(r.headers.get("Content-Length") or 0)
                while True:
                    b = r.read(1 << 18)
                    if not b: break
                    f.write(b); h.update(b); size += len(b)
            if expect and size != expect: raise IOError("short read")
            break
        except Exception as e:
            print(qt, "retry", attempt, e, flush=True); time.sleep(15)
    else:
        print(qt, "DOWNLOAD FAILED", flush=True); return
    same = h.hexdigest() == MAN[qt]["zip_sha256"]
    with plock:
        try:
            with zipfile.ZipFile(zp) as z:
                m = [n for n in z.namelist() if re.search(r"(^|/)DEMO\d\dq\d(_new)?\.txt$", n, re.I)][0]
                with z.open(m) as f:
                    df = pd.read_csv(f, sep="$", dtype=str, quoting=3, encoding="latin-1", on_bad_lines="skip", index_col=False, engine="c", low_memory=False, keep_default_na=False)
            df.columns = [c.strip().lower() for c in df.columns]; o = pd.DataFrame(index=df.index)
            for name, alts in WANT.items():
                src = next((a for a in alts if a in df.columns), None)
                col = df[src].str.strip().str.upper() if src else pd.Series([""] * len(df), index=df.index)
                o[name] = pd.to_numeric(col, errors="coerce").astype("Int64") if name in ("pid", "caseid") else col
            o = o[o.pid.notna()]; o.to_parquet(out + ".part", compression="zstd", index=False); os.replace(out + ".part", out)
        except Exception as e:
            print(qt, "PARSE FAILED", repr(e)[:200], flush=True); return
        finally:
            if os.path.exists(zp): os.remove(zp)
    with lock:
        done[qt] = dict(quarter=qt, rows=len(o), zip_sha256=h.hexdigest(), same_file_as_first_pass=int(same))
        with open(LOG, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["quarter", "rows", "zip_sha256", "same_file_as_first_pass"], delimiter="\t"); w.writeheader()
            for k in sorted(done): w.writerow(done[k])
        print(qt, len(o), "same file" if same else "DIFFERENT FILE", f"[{len(done)}]", flush=True)
shutil.rmtree(TMP, ignore_errors=True); os.makedirs(TMP)
with ThreadPoolExecutor(14) as ex: list(ex.map(work, sorted(MAN, reverse=True)))
print("finished", len(done), "of", len(MAN))
