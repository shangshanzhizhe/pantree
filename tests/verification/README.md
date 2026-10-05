# Randomized verification

Scripts used to check the fixes on this branch. pytest does not collect them; they are meant for reviewers and could be
dropped before merging.

| Script | What it does |
|---|---|
| `validate_simplified.py IN.gfa OUT.gfa` | Independent check of a simplified GFA. Each output segment, read as its `oi` list, must spell a path of the input graph. Each output link must join the boundary segments of its ends through an input link. No input segment may appear twice. |
| `fuzz3.py N OUTDIR` | Generates small random graphs, most with cyclic reference walks. Each graph goes through `from_gfa`, `simplify_subgraph`, `write_GFA`, `write_vcf` and `consolidate` with both DFS methods, checking the tree, strand and path invariants. The smallest failing graph of each kind is written to `OUTDIR`. |
| `dump2.py N OUTDIR` + `compare2.py OLD NEW` | Random graphs with an acyclic reference walk and no self edges, run through both DFS methods. Compares trees, VCFs (with inversion RC/AN/CR normalised to `.`) and simplified topology between two versions, and validates both simplified outputs. |
| `fuzz.py`, `fuzz_dump.py` | Graph generators used by the scripts above. |

Run the scripts from this directory. The version under test is whichever `pantree` comes first on `PYTHONPATH`; avoid
running from the repository root, where the local `pantree/` always wins.

```bash
cd tests/verification
export PYTHONHASHSEED=0                        # fixed seed so that both versions see identical inputs

# this branch: all invariants on random graphs, without and with self edges
PYTHONPATH=../..:. python fuzz3.py 2000 out_noself
ALLOW_SELF_EDGES=1 PYTHONPATH=../..:. python fuzz3.py 2000 out_self

# original vs this branch on acyclic reference walks
git -C ../.. worktree add /tmp/pantree-main main
PYTHONPATH=/tmp/pantree-main:. python dump2.py 3000 d_main
PYTHONPATH=../..:. python dump2.py 3000 d_branch
python compare2.py d_main d_branch
```

Expected:
- `fuzz3.py` reports only `ok`.
- `compare2.py` reports identical trees, VCFs and simplified topology, and `PATCHED has invalid segments` is 0.
- The original typically has invalid segments in about two thirds of the graphs.
