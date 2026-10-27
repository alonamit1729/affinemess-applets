#!/bin/bash
cd /tmp/claude-0/-home-claude/27d8f328-c517-545a-96bf-782e414fcf24/scratchpad/symtest
for i in $(seq 9 112); do timeout 600 python3 msolve_check.py 6,7 $i 2 || echo "6,7 $i TIMEOUT/FAIL"; done
