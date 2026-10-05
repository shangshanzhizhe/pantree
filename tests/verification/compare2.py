"""Compare dump2.py outputs of the original and patched code, allowing only the intended differences."""
import json
import os
import sys
from collections import Counter

from validate_simplified import validate, read_gfa

FLIP = {'+': '-', '-': '+'}
old, new = sys.argv[1], sys.argv[2]


def norm_vcf(path):
    """Drop headers; for inversions report RC, AN and CR as '.', which is what the patch now writes."""
    rows = []
    for line in open(path):
        if line.startswith('##'):
            continue
        p = line.rstrip('\n').split('\t')
        if not line.startswith('#') and 'VT=INV' in p[7]:
            info = [f'{k}=.' if k in ('RC', 'AN') else f'{k}={v}'
                    for k, _, v in (kv.partition('=') for kv in p[7].split(';'))]
            p[7] = ';'.join(info)
            for i in range(9, len(p)):
                gt, cr, ca = p[i].split(':')
                p[i] = ':'.join([gt, ','.join('.' for _ in cr.split(',')), ca])
        rows.append(p)
    return rows


def canonical_topology(path):
    segs, links, oi = read_gfa(path)
    key, flip = {}, {}
    for sid in segs:
        ids = tuple(oi.get(sid, []))
        if not ids:  # artificial terminus segment without provenance
            continue
        rev = tuple(reversed(ids))
        key[sid], flip[sid] = min(ids, rev), rev < ids
    canon_links = set()
    for a, oa, b, ob in links:
        if a not in key or b not in key:
            continue
        a2, oa2 = key[a], (FLIP[oa] if flip[a] else oa)
        b2, ob2 = key[b], (FLIP[ob] if flip[b] else ob)
        edge = (a2, oa2, b2, ob2)
        canon_links.add(min(edge, (b2, FLIP[ob2], a2, FLIP[oa2])))
    return set(key.values()), canon_links


stats = Counter()
examples = {}
cases = sorted({f.split('.')[0] for f in os.listdir(old) if f.endswith('.gfa') and f.count('.') == 1})
for c in cases:
    gfa = os.path.join(old, c + '.gfa')
    for method in ('max_weight', 'contiguous'):
        o, n = os.path.join(old, f'{c}.{method}'), os.path.join(new, f'{c}.{method}')
        of, nf = os.path.exists(o + '.FAILED'), os.path.exists(n + '.FAILED')
        if of or nf:
            k = f'{method}: ' + ('fail in both' if of and nf else 'fixed' if of else 'NEWLY BROKEN')
            if of:
                k += ' (' + open(o + '.FAILED').readline().strip() + ')'
            stats[k] += 1
            examples.setdefault(k, c)
            continue
        stats[f'{method}: tree ' + ('same' if open(o + '.tree').read() == open(n + '.tree').read() else 'DIFFERENT')] += 1
        stats[f'{method}: vcf ' + ('same' if norm_vcf(o + '.vcf') == norm_vcf(n + '.vcf') else 'DIFFERENT')] += 1
        same_topology = canonical_topology(o + '.simplified.gfa') == canonical_topology(n + '.simplified.gfa')
        stats[f'{method}: simplified topology ' + ('same' if same_topology else 'DIFFERENT')] += 1
        vo, vn = validate(gfa, o + '.simplified.gfa'), validate(gfa, n + '.simplified.gfa')
        stats[f'{method}: original has invalid segments'] += bool(vo['invalid_segments'])
        stats[f'{method}: PATCHED has invalid segments'] += bool(vn['invalid_segments'] or vn['invalid_links'])
for k, v in sorted(stats.items()):
    print(f'  {v:6d}  {k}' + (f'   (e.g. {examples[k]})' if k in examples else ''))
