"""Random small bidirected graphs with cyclic reference walks; run pantree end to end and classify failures."""
import os
import random
import sys
import tempfile
import traceback
import warnings
from collections import Counter

warnings.simplefilter("ignore")
import logging
logging.disable(logging.CRITICAL)
from pantree.graph import PangenomeGraph

FLIP = {'>': '<', '<': '>'}


def random_walk(rng, n_segments, length, ref=None):
    if ref is not None and rng.random() < 0.6:
        # haplotype derived from the reference: drop some steps / revisits
        walk = [s for s in ref if rng.random() < 0.75]
        return walk or ref[:1]
    walk = [(rng.choice('>' * 5 + '<'), rng.randint(1, n_segments))]
    while len(walk) < length:
        walk.append((rng.choice('>' * 4 + '<'), rng.randint(1, n_segments)))
    return walk


def make_case(seed, max_segments=10):
    rng = random.Random(seed)
    n = rng.randint(3, max_segments)
    # reference: mostly new segments in order, with revisits (same or opposite orientation)
    ref, next_new = [], 1
    for _ in range(rng.randint(4, 2 * n + 2)):
        if next_new <= n and (rng.random() < 0.6 or not ref):
            ref.append(('>' if rng.random() < 0.9 else '<', next_new))
            next_new += 1
        else:
            o, s = rng.choice(ref)
            ref.append((o if rng.random() < 0.5 else FLIP[o], s))
    walks = {'REF': ref}
    for h in range(rng.randint(1, 3)):
        walks[f'H{h}'] = random_walk(rng, n, rng.randint(2, 2 * n), ref)
    used = sorted({s for w in walks.values() for _, s in w})
    seqs = {s: ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 6))) for s in used}
    links = set()
    for w in walks.values():
        for (o1, a), (o2, b) in zip(w, w[1:]):
            e = (a, o1, b, o2)
            c = (b, FLIP[o2], a, FLIP[o1])
            if e not in links and c not in links:
                links.add(e)
    lines = ['H\tVN:Z:1.1']
    lines += [f'S\t{s}\t{seqs[s]}' for s in used]
    sign = {'>': '+', '<': '-'}
    lines += [f'L\t{a}\t{sign[o1]}\t{b}\t{sign[o2]}\t0M' for a, o1, b, o2 in sorted(links)]
    for name, w in walks.items():
        lines.append(f'W\t{name}\t0\tchr1\t0\t0\t' + ''.join(f'{o}{s}' for o, s in w))
    return '\n'.join(lines) + '\n', walks


def no_self_edges(walks):
    return all(a != b for w in walks.values() for (_, a), (_, b) in zip(w, w[1:]))


def run_case(text, min_len, tmpdir):
    path = os.path.join(tmpdir, 'g.gfa')
    with open(path, 'w') as f:
        f.write(text)
    G = PangenomeGraph.from_gfa(path, ref_name='REF')
    G.simplify_subgraph(minimum_allele_length=min_len)
    G.write_vcf(path, os.path.join(tmpdir, 'out.vcf'), 'chr1')


def classify(exc):
    tb = traceback.extract_tb(exc.__traceback__)
    frames = [f for f in tb if '/pantree/' in f.filename and 'site-packages' not in f.filename]
    where = f'{os.path.basename(frames[-1].filename)}:{frames[-1].name}' if frames else '?'
    return f'{type(exc).__name__} @ {where}'


def main(n_cases, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    counts, smallest = Counter(), {}
    cyclic = 0
    with tempfile.TemporaryDirectory() as tmp:
        seed = -1
        while sum(counts.values()) < n_cases:
            seed += 1
            text, walks = make_case(seed)
            if not ALLOW_SELF_EDGES and not no_self_edges(walks):
                continue
            segs = [s for _, s in walks['REF']]
            cyclic += len(segs) != len(set(segs))
            try:
                run_case(text, random.Random(seed).randint(1, 30), tmp)
                counts['ok'] += 1
            except Exception as e:  # noqa: BLE001
                key = classify(e)
                counts[key] += 1
                if key not in smallest or len(text) < len(smallest[key][1]):
                    smallest[key] = (seed, text)
    print(f'cases: {n_cases}, cyclic reference: {cyclic}')
    for k, v in counts.most_common():
        print(f'  {v:6d}  {k}')
    for k, (seed, text) in smallest.items():
        fn = os.path.join(out_dir, k.split(' @ ')[0] + '.' + k.split(':')[-1] + '.gfa')
        with open(fn, 'w') as f:
            f.write(text)
        print(f'smallest for {k}: seed {seed} -> {fn}')


ALLOW_SELF_EDGES = os.environ.get('ALLOW_SELF_EDGES') == '1'

if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2])
