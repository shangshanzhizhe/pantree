"""Write simplify/VCF outputs for random graphs whose reference walk is acyclic, for comparing two pantree versions."""
import os
import random
import sys
import traceback

import fuzz
from pantree.graph import PangenomeGraph
from pantree.write_gfa import write_GFA


def make_acyclic_case(seed):
    text, walks = fuzz.make_case(seed)
    rng = random.Random(seed + 10**6)
    ref, seen = [], set()
    for o, s in walks['REF']:  # drop revisits from the reference only
        if s not in seen:
            ref.append((o, s))
            seen.add(s)
    walks['REF'] = ref
    lines = [l for l in text.splitlines() if not l.startswith(('L', 'W'))]
    links = set()
    for w in walks.values():
        for (o1, a), (o2, b) in zip(w, w[1:]):
            e, c = (a, o1, b, o2), (b, fuzz.FLIP[o2], a, fuzz.FLIP[o1])
            if e not in links and c not in links:
                links.add(e)
    sign = {'>': '+', '<': '-'}
    lines += [f'L\t{a}\t{sign[o1]}\t{b}\t{sign[o2]}\t0M' for a, o1, b, o2 in sorted(links)]
    lines += [f'W\t{n}\t0\tchr1\t0\t0\t' + ''.join(f'{o}{s}' for o, s in w) for n, w in walks.items()]
    return '\n'.join(lines) + '\n', walks, rng.randint(1, 30)


def main(n_cases, out_dir, allow_self):
    os.makedirs(out_dir, exist_ok=True)
    done, seed = 0, -1
    while done < n_cases:
        seed += 1
        text, walks, min_len = make_acyclic_case(seed)
        if not allow_self and not fuzz.no_self_edges(walks):
            continue
        done += 1
        gfa = os.path.join(out_dir, f'{seed}.gfa')
        with open(gfa, 'w') as f:
            f.write(text)
        try:
            G = PangenomeGraph.from_gfa(gfa, ref_name='REF')
            write_GFA(G.simplify_subgraph(minimum_allele_length=min_len), os.path.join(out_dir, f'{seed}.simplified.gfa'))
            G.write_vcf(gfa, os.path.join(out_dir, f'{seed}.vcf'), 'chr1')
        except Exception as e:  # noqa: BLE001
            with open(os.path.join(out_dir, f'{seed}.FAILED'), 'w') as f:
                f.write(fuzz.classify(e) + '\n' + traceback.format_exc())


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], os.environ.get('ALLOW_SELF_EDGES') == '1')
