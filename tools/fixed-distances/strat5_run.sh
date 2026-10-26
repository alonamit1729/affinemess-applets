#!/bin/bash
cd /tmp/claude-0/-home-claude/27d8f328-c517-545a-96bf-782e414fcf24/scratchpad
python3 strat.py 5 >> strat5.log 2>&1
python3 strat_analyze.py 5 > strat5_analysis.txt 2>&1
