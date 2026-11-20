"""Write cells-table.tex (the longtable of solid cells) for table-methods.tex from ../data/applet_data.json.
Run from this folder:  python3 make_tables.py"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, '..', 'data', 'applet_data.json')))
WHERE = {   # where the completeness proof comes from (section of table-methods.tex, or the literature)
 '3,1': 'trivial', '3,2': 'trivial', '3,3': 'trivial',
 '4,2': r'\cite{ES}; \S\ref{sec:small}', '4,3': r'\S\ref{sec:small}', '4,4': r'\S\ref{sec:small}', '4,5': r'\S\ref{sec:small}',
 '4,6': r'\S\ref{sec:small}', '5,2': r'\cite{ES}; \S\ref{sec:small}', '5,3': r'\cite{Sh04}; \S\ref{sec:small}',
 '5,4': r'\S\ref{sec:small}', '5,5': r'\S\ref{sec:small}', '5,6': r'\S\ref{sec:small}', '5,7': r'\S\ref{sec:small}',
 '6,3': r'\cite{Sh04}; \S\ref{sec:published}', '7,3': r'\cite{EF,Sh04}; \S\ref{sec:published}',
 '7,4': r'\cite{LW}; \S\ref{sec:published}', '8,4': r'\cite{Sh08}; \S\ref{sec:published}', '9,4': r'\cite{EF}; \S\ref{sec:published}',
 '11,5': r'\cite{W11}; \S\ref{sec:published}', '12,5': r'\cite{EF,Sh08}; \S\ref{sec:published}',
 '13,6': r'\cite{W12,SO}; \S\ref{sec:published}', '11,6': r'\cite{AA}; \S\ref{sec:ten}',
 '10,5': r'\S\ref{sec:ten}; \cite[Cor.~7]{AC}', '10,6': r'\S\ref{sec:ten}; \cite[Cor.~7]{AC}', '12,6': r'\S\ref{sec:ten}',
 '14,7': r'\cite{AC14}', '15,7': r'\cite[Thm.~16]{AC}', '16,7': r'\cite[Thm.~1]{AC}',
}
for k in ('6,4', '6,5', '6,6', '7,5', '7,6', '8,5', '8,6', '9,5', '9,6'): WHERE[k] = r'\S\ref{sec:middle}'
rows = []
for k in sorted(D['cells'], key=lambda s: tuple(map(int, s.split(',')))):
    c = D['cells'][k]
    if c['status'] != 'exact': continue
    dims = {int(a): b['lo'] for a, b in c['dims'].items()}
    rig = dims.get(0, 0)
    fam = ', '.join(f'{v}\\,({a})' for a, v in sorted(dims.items()) if a > 0) or '--'
    mark = {'†': r'$\dagger$', '‡': r'$\ddagger$'}.get(''.join(c['marks']), '')
    n, d = k.split(',')
    rows.append(f'$({n},{d})$ & {rig} & {fam} & {mark} & {WHERE[k]} \\\\')
HEAD = (r'\begin{longtable}{@{}lrllp{0.42\textwidth}@{}}' '\n'
        r'\caption{The solid cells.}\label{tab:cells}\\' '\n'
        r'\toprule cell & rigid & families & mark & proof of completeness \\ \midrule' '\n' r'\endfirsthead' '\n'
        r'\toprule cell & rigid & families & mark & proof of completeness \\ \midrule' '\n' r'\endhead' '\n')
open(os.path.join(HERE, 'cells-table.tex'), 'w').write(HEAD + '\n'.join(rows) + '\n' + r'\bottomrule' '\n' r'\end{longtable}' '\n')
print(len(rows), 'rows')
