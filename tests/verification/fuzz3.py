"""End-to-end fuzzing of every pantree command with output invariants (run with the dev copy on sys.path)."""
import os
import random
import sys
import tempfile
import warnings
from collections import Counter

warnings.simplefilter("ignore")
import logging
logging.disable(logging.CRITICAL)

import fuzz
from validate_simplified import validate
from pantree.graph import PangenomeGraph
from pantree.dfs import dfs_methods
from pantree.utils import node_complement, sequence_complement, edge_complement
from pantree.write_gfa import write_GFA
from pantree.walk_variants import process_haplotype_variants


class Violation(Exception):
    pass


def check_tree(G):
    T = G.reference_tree
    roots = [n for n in T if T.in_degree(n) == 0]
    if roots != ['+_terminus_+']:
        raise Violation(f'tree roots {roots[:3]}')
    if T.number_of_nodes() and max(d for _, d in T.in_degree()) > 1:
        raise Violation('tree node with two parents')
    missing = [n for n, d in G.nodes(data='direction') if d == 1 and n not in T]
    if missing:
        raise Violation('positive node missing from tree')


def check_simplified(S):
    for n in S:
        if S.nodes[n].get('sequence', '') != sequence_complement(S.nodes[node_complement(n)].get('sequence', '')):
            raise Violation('strands disagree after simplify')
    for e in S.edges:
        if not S.has_edge(*edge_complement(e)):
            raise Violation('edge without complement after simplify')


RENAME = {'REF': 'refs', 'H0': 's0', 'H1': 's1', 'H2': 's2'}  # sample names must not clash with VCF columns


def rename_samples(text):
    lines = []
    for line in text.splitlines():
        p = line.split('\t')
        if p[0] == 'W':
            p[1] = RENAME[p[1]]
        lines.append('\t'.join(p))
    return '\n'.join(lines) + '\n'


def run(text, walks, min_len, method, tmp):
    text = rename_samples(text)
    gfa = os.path.join(tmp, 'g.gfa')
    with open(gfa, 'w') as f:
        f.write(text)
    priority = {'refs': 0, 's0': 1} if method == 'contiguous' else None
    G = PangenomeGraph.from_gfa(gfa, ref_name='refs', dfs_method=dfs_methods[method], priority_dict=priority)
    check_tree(G)
    S = G.simplify_subgraph(minimum_allele_length=min_len)
    check_simplified(S)
    out = os.path.join(tmp, 's.gfa')
    write_GFA(S, out)
    v = validate(gfa, out)
    if v['invalid_segments'] or v['invalid_links'] or v['duplicated_input_segments']:
        raise Violation(f"simplified GFA invalid: {v['invalid_segments']} segments, {v['invalid_links']} links")
    vcf = os.path.join(tmp, 'o.vcf')
    G.write_vcf(gfa, vcf, 'chr1')
    if 'H0' in walks:
        process_haplotype_variants(vcf, 's0', 0, os.path.join(tmp, 'c.vcf'))


def main(n_cases, out_dir, allow_self, cyclic_only):
    os.makedirs(out_dir, exist_ok=True)
    counts, smallest = Counter(), {}
    with tempfile.TemporaryDirectory() as tmp:
        seed, done = -1, 0
        while done < n_cases:
            seed += 1
            text, walks = fuzz.make_case(seed)
            if not allow_self and not fuzz.no_self_edges(walks):
                continue
            segs = [s for _, s in walks['REF']]
            if cyclic_only and len(segs) == len(set(segs)):
                continue
            done += 1
            min_len = random.Random(seed).randint(1, 30)
            for method in ('max_weight', 'contiguous'):
                try:
                    run(text, walks, min_len, method, tmp)
                    counts[f'{method}: ok'] += 1
                except Violation as e:
                    key = f'{method}: VIOLATION {e}'
                except Exception as e:  # noqa: BLE001
                    key = f'{method}: {fuzz.classify(e)}'
                else:
                    continue
                counts[key] += 1
                if key not in smallest or len(text) < len(smallest[key][1]):
                    smallest[key] = (seed, text, min_len)
    for k, v in sorted(counts.items()):
        print(f'  {v:6d}  {k}')
    for i, (k, (seed, text, min_len)) in enumerate(sorted(smallest.items())):
        fn = os.path.join(out_dir, f'case{i}.seed{seed}.min{min_len}.gfa')
        with open(fn, 'w') as f:
            f.write(text)
        print(f'  smallest for [{k}]: {fn}')


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], os.environ.get('ALLOW_SELF_EDGES') == '1',
         os.environ.get('CYCLIC_ONLY') == '1')
