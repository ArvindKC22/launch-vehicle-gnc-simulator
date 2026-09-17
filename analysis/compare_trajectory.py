#!/usr/bin/env python3
"""Compare simulator summary metrics with an OpenRocket/RASAero CSV export."""
import argparse, csv, math

def read(path):
    with open(path, newline="") as f:
        return next(csv.DictReader(f))

def value(row, names):
    for name in names:
        if name in row and row[name] not in ("", "nan", "NaN"):
            return float(row[name])
    raise KeyError("missing one of: " + ", ".join(names))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("simulator_csv")
    ap.add_argument("reference_csv", help="single-row CSV exported from OpenRocket or RASAero")
    a=ap.parse_args(); sim=read(a.simulator_csv); ref=read(a.reference_csv)
    metrics=[("apogee_m",["apogee_m","apogee","Apogee (m)"]),
             ("range_m",["impact_x_m","range_m","range","Max. altitude range (m)"])]
    print("metric,simulator,reference,absolute_error,percent_error")
    for label,names in metrics:
        try: s=value(sim,names); r=value(ref,names)
        except KeyError: continue
        err=s-r; pct=100.0*err/r if r else math.nan
        print(f"{label},{s:.9g},{r:.9g},{err:.9g},{pct:.6g}")
if __name__ == "__main__": main()
