"""Independent check of a simplified GFA against its input GFA.

Each output segment, read as its `oi` list of input segments, must be a path of the input graph whose oriented
sequences concatenate to the output sequence; each output link must join the boundary segments of its ends by an
input link; and no input segment may appear in two output segments.
"""
import json
import sys

COMP = str.maketrans('ACGTRYKMBVDHNacgtrykmbvdhn', 'TGCAYRMKVBHDNtgcayrmkvbhdn')
FLIP = {'+': '-', '-': '+'}


def rc(s):
    return s.translate(COMP)[::-1]


def read_gfa(path):
    segs, links, oi = {}, set(), {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if p[0] == 'S':
                segs[p[1]] = p[2] if p[2] != '*' else ''
                for tag in p[3:]:
                    if tag.startswith('oi:J:'):
                        oi[p[1]] = [str(x) for x in json.loads(tag[5:])]
            elif p[0] == 'L':
                a, oa, b, ob = p[1:5]
                links.add((a, oa, b, ob))
                links.add((b, FLIP[ob], a, FLIP[oa]))
    return segs, links, oi


def oriented(seq, o):
    return seq if o == '+' else rc(seq)


def endpoint_options(pieces, seq, in_segs, in_links):
    """Feasible (first oriented piece, last oriented piece) pairs, or empty if the segment is not a valid path."""
    if not pieces:
        return {(None, None)} if seq == '' else set()
    states = {}  # orientation of current piece -> set of first oriented pieces
    pos = 0
    for i, pid in enumerate(pieces):
        s = in_segs[pid]
        new = {}
        for o in '+-':
            if seq[pos:pos + len(s)] != oriented(s, o):
                continue
            if i == 0:
                new.setdefault(o, set()).add((pid, o))
            else:
                for po, firsts in states.items():
                    if (pieces[i - 1], po, pid, o) in in_links:
                        new.setdefault(o, set()).update(firsts)
        states, pos = new, pos + len(s)
        if not states:
            return set()
    if pos != len(seq):
        return set()
    return {(first, (pieces[-1], o)) for o, firsts in states.items() for first in firsts}


def validate(in_path, out_path):
    in_segs, in_links, _ = read_gfa(in_path)
    out_segs, out_links, out_oi = read_gfa(out_path)
    options, bad_segments, seen, duplicated = {}, [], {}, []
    for sid, seq in out_segs.items():
        pieces = out_oi.get(sid, [])
        for pid in pieces:
            if pid in seen:
                duplicated.append(pid)
            seen[pid] = sid
        options[sid] = endpoint_options(pieces, seq, in_segs, in_links)
        if not options[sid]:
            bad_segments.append(sid)
    bad_links = 0
    for a, oa, b, ob in out_links:
        if not options.get(a) or not options.get(b) or not out_oi.get(a) or not out_oi.get(b):
            continue  # invalid segment already counted, or an artificial terminus segment
        ok = False
        for first_a, last_a in options[a]:
            end_a = last_a if oa == '+' else (first_a[0], FLIP[first_a[1]])
            for first_b, last_b in options[b]:
                start_b = first_b if ob == '+' else (last_b[0], FLIP[last_b[1]])
                if (end_a[0], end_a[1], start_b[0], start_b[1]) in in_links:
                    ok = True
        bad_links += not ok
    return {'segments': len(out_segs), 'invalid_segments': len(bad_segments), 'invalid_links': bad_links // 1,
            'duplicated_input_segments': len(duplicated), 'examples': bad_segments[:3]}


if __name__ == '__main__':
    print(validate(sys.argv[1], sys.argv[2]))
