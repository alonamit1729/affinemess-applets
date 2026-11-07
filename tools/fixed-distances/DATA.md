# Large data files (not in the published tree)

GitHub Pages publishes everything on `main`, so files that only matter for *reproducing* the work are kept out of
the tree (2026-11-06). Nothing was lost: each file is still in git history, byte for byte, at the pinned commit
`ba9fb77dc1f62e81d4d550a6bc5405a674aa97c3`, and `fetch_data.py` restores it and checks its SHA-256:

    python3 tools/fixed-distances/fetch_data.py            # all eleven (about 84 MB)
    python3 tools/fixed-distances/fetch_data.py cells/rigid_6_6.json
    python3 tools/fixed-distances/fetch_data.py --list

The restored files are git-ignored. Rule used: files over 1 MB that no page fetches. Browse them at
<https://github.com/alonamit1729/affinemess-applets/tree/ba9fb77dc1f62e81d4d550a6bc5405a674aa97c3/tools/fixed-distances>.

| file | size | SHA-256 | made by | used by |
|---|---|---|---|---|
| `data/exact_cells2.json` | 22.9 MB | `7d70a04236846f13b4efd75fdf17986f4bf28ea05f2a1b57f56123f6dc597251` | `cells/export_cells.py` (from `cells/rigid_*`, `cells/fam_*`) | `build_data2.py` (exact cells of the table and catalogs) |
| `data/strat_4.json` | 14.3 MB | `0d37e81029af6ad487622778e5277553be113ff0a99498aae8bf1be15fe723b4` | `strat.py 4` | `strat_analyze.py` (writes the small `strat_graph_*.json`) |
| `data/strat_5.json` | 8.9 MB | `f995d94fa944e90884c0ab08b950f2153800aac867739e7b32458a2c5d508caa` | `strat.py 5` | `strat_analyze.py` |
| `data/degen_5.json` | 2.6 MB | `caa05ef4ff65e3f188781f6b36cade73457f3d806c17f1962df9b3d568fdb8b0` | `degen.py 5` | `degen_analyze2.py` (writes `degen_poset_5.json`) |
| `data/degen_s_5.json` | 1.8 MB | `e851454a3ebec9d6228aebd9f8aac850b5db6850df619705a19ca663334e1a3d` | `degen_sample.py 5` | `degen_analyze2.py` |
| `data/degen_s_6.json` | 1.7 MB | `a22394a29112567bf5716818e54a092ee88e183573d022423bef0586a10fc5f3` | `degen_sample.py 6` | `degen_analyze2.py` |
| `data/sat56.json` | 1.5 MB | `dbd653cb2a39ef7f8a99964a251bba8941dcd25c527494be8a84c331dbb53d96` | `sat_n.py 5 6` | `strat_n.py`, `degen_n.py` |
| `cells/cert_6_6.json` | 13.6 MB | `45cb1ebc5c176e4aa55916d497a77299022e7c8212757cf4a93777a5df671491` | `cells/cellcert.py 6 6` (about 6 h) | `cells/summary.py`, `cells/export_cells.py` |
| `cells/rigid_6_6.json` | 11.7 MB | `1d0b51ab6ffa6d726f3aad43ea8214834f0b06cda5f2f13a1433f60089bda1e2` | `cells/summary.py 6 6` | `cells/export_cells.py` |
| `cells/cert_7_6.json` | 2.7 MB | `9db7c3ab6eb17d66fae9a68c877b845fc522767e78e2053b9da42f6351848c42` | `cells/cellcert.py 7 6` | `cells/summary.py`, `cells/export_cells.py` |
| `cells/rigid_7_6.json` | 2.2 MB | `40a426d98f94107494f6e198cd7e297e85312225af0e5aa27a08846af1fee0ac` | `cells/summary.py 7 6` | `cells/export_cells.py` (maximality of the (6,6) sets) |

## What stays in the tree, and why

- `data/applet_data.json` (3.6 MB, *slim* since 2026-11-07): the master data of the applet. The items of the four cells with
  more than 1000 configurations are not repeated in it, only their counts and a `remote` pointer; they live in
  `data/catalogs/`. `datautil.load('data/applet_data.json')` returns the full structure (it was 47 MB), and the
  scripts that need every item (`cells/export_cells.py`, `symtest/`, `strat_analyze.py`, `g7fin/catalog_match.py`, Distance
  Hunt's `canon.py`, ...) call it. `build_html.py` writes the catalogs and slims the file; it is idempotent and rewrites a
  file only when its bytes change.
- `data/catalogs/*.json` (38 MB; the cards keep their numbers, the page shows families first): fetched on demand by the live page when a cell with more than 1000 configurations
  is opened. `catalogs/6_6.json` alone holds the 44 242 + 757 items of the (6,6) cell.
- `cells/fam_*.json` and the smaller `cells/cert_*`, `cells/rigid_*` (under 1 MB each): the classification results
  of the other cells.

The (6,6) classification itself is therefore complete in the tree: `cells/fam_6_6.json` (the 757 families) and, for
the 44 242 rigid configurations with coordinates and labels, `data/catalogs/6_6.json`; `cells/rigid_6_6.json` and
`cells/cert_6_6.json` are the same results in the form `summary.py` and `export_cells.py` read, plus the per-pattern
certificate.

`build_data2.py` needs `data/exact_cells2.json`. It can be fetched, or rebuilt from the `cells/` results with
`python3 cells/export_cells.py n d [n d ...]` for the cells that have a `cells/rigid_n_d.json` (the (6,6) cell takes its maximality labels from `cells/rigid_7_6.json`; fetch it first).
