#!/bin/bash
# sequential, single-process pipeline (temp files are pid-unique; nothing else may use msolve concurrently)
cd /tmp/claude-0/-home-claude/27d8f328-c517-545a-96bf-782e414fcf24/scratchpad
set -x
python3 certify5d.py 5 >> cert55_redo.log 2>&1
python3 refold_noreal.py >> hyg.log 2>&1
for d in 5 6; do for i in $(python3 -c "import json;print(' '.join(str(r['idx']) for r in json.load(open('cert5$d.json')) if r.get('unresolved')))"); do [ -f stubborn_${d}_$i.json ] || python3 stubborn.py $d $i 2400 > stubborn_${d}_$i.log 2>&1; done; done
for d in 7 6 5; do python3 hidden.py $d >> hidden$d.log 2>&1; done
for d in 7 6 5; do python3 hidden2.py $d 1200 >> hidden2_$d.log 2>&1; done
for d in 7 6 5; do python3 hidden2b.py $d > hidden2b_$d.log 2>&1; python3 hidden2c.py $d > hidden2c_$d.log 2>&1; python3 hidden3.py $d > hidden3_$d.log 2>&1; done
echo alldone > hygiene.done
