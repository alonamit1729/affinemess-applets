#!/usr/bin/env python3
"""Export a level's pattern list to the g7fin/patterns_{n}.json format
   [{"name": ..., "k": ..., "col": {"i,j": colour}, "sources": [...]}]
Usage: python3 export.py n [--levels DIR] [--which L|raw] [--out FILE]
"""
import argparse, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from patlib import col_to_dict

ap = argparse.ArgumentParser()
ap.add_argument('n', type=int)
ap.add_argument('--levels', default=os.path.join(HERE, 'levels'))
ap.add_argument('--which', default='L', choices=['L', 'raw'])
ap.add_argument('--out', default=None)
args = ap.parse_args()
n = args.n
pats = json.load(open(os.path.join(args.levels, f'{args.which}{n}.json')))['patterns']
out = [{'name': f'gen{n}_{i}', 'k': max(c) + 1, 'col': col_to_dict(n, c), 'sources': [f'gen/{args.which}{n}']}
       for i, c in enumerate(pats)]
fn = args.out or os.path.join(args.levels, f'patterns_{n}_{args.which}.json')
json.dump(out, open(fn, 'w'))
print(f'{len(out)} patterns written to {fn}')
