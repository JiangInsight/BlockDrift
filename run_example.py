import csv
from itertools import product
from pathlib import Path
import struct
import subprocess
import tempfile

root = Path(__file__).resolve().parent
with (root/"data/facts.csv").open(newline="") as stream:
    facts = [tuple(map(int,r)) for r in list(csv.reader(stream))[1:]]
bins, keys = max(t for t,g in facts)+1, max(g for t,g in facts)+1
windows = [(a,b) for a in range(bins+1) for b in range(a,bins+1)]
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    packed, schedule, oracle = (work/n for n in ("facts.bin","schedule.csv","oracle.csv"))
    packed.write_bytes(b"".join(struct.pack("<ii",t,g) for t,g in facts))
    with schedule.open("w",newline="") as requests, oracle.open("w",newline="") as truth:
        q, r = csv.writer(requests), csv.writer(truth)
        q.writerow(["round","step","qid","a","b","c","d","k"])
        r.writerow(["qid","rank","key","delta"])
        for i, ((a,b),(c,d),k) in enumerate(product(windows,windows,[1,3,6])):
            q.writerow([0,i,i,a,b,c,d,k])
            delta = [0]*keys
            for time,key in facts:
                delta[key] += int(c<=time<d) - int(a<=time<b)
            ranked = sorted((g for g in range(keys) if delta[g]),key=lambda g:(-abs(delta[g]),g))
            for rank,g in enumerate(ranked[:k],1):
                r.writerow([i,rank,g,delta[g]])
    for method in ("split","drift","prefixfull","prefix256"):
        subprocess.run([root/"build/blockdrift.exe",packed,str(bins),str(keys),method,
                        schedule,oracle,work/"observations.csv",work/"build.csv"],check=True)
print(f"Four execution paths agree with direct counts on {3*len(windows)**2} ranked requests each.")
