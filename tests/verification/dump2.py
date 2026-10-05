"""Dump trees, simplified GFAs and VCFs for random graphs with an acyclic reference walk (no self edges)."""
import os
import sys
import traceback
import warnings

warnings.simplefilter("ignore")
import logging
logging.disable(logging.CRITICAL)

import fuzz
from fuzz_dump import make_acyclic_case
from fuzz3 import rename_samples
from pantree.graph import PangenomeGraph
from pantree.dfs import dfs_methods
from pantree.write_gfa import write_GFA


def main(n_cases, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    done, seed = 0, -1
    while done < n_cases:
        seed += 1
        text, walks, min_len = make_acyclic_case(seed)
        if not fuzz.no_self_edges(walks):
            continue
        done += 1
        text = rename_samples(text)
        gfa = os.path.join(out_dir, f'{seed}.gfa')
        with open(gfa, 'w') as f:
            f.write(text)
        for method in ('max_weight', 'contiguous'):
            stem = os.path.join(out_dir, f'{seed}.{method}')
            try:
                priority = {'refs': 0, 's0': 1} if method == 'contiguous' else None
                G = PangenomeGraph.from_gfa(gfa, ref_name='refs', dfs_method=dfs_methods[method], priority_dict=priority)
                with open(stem + '.tree', 'w') as f:
                    f.write('\n'.join(sorted(f'{u}\t{v}' for u, v in G.reference_tree.edges())) + '\n')
                write_GFA(G.simplify_subgraph(minimum_allele_length=min_len), stem + '.simplified.gfa')
                G.write_vcf(gfa, stem + '.vcf', 'chr1')
            except Exception as e:  # noqa: BLE001
                with open(stem + '.FAILED', 'w') as f:
                    f.write(fuzz.classify(e) + '\n' + traceback.format_exc())


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2])
