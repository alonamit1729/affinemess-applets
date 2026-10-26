#!/bin/bash
cd /tmp/claude-0/-home-claude/27d8f328-c517-545a-96bf-782e414fcf24/scratchpad
python3 degen_sample.py 6 >> degen_s6.log 2>&1
python3 degen_sample.py 7 >> degen_s7.log 2>&1
python3 degen_sample.py 5 >> degen_s5.log 2>&1
