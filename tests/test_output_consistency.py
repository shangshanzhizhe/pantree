"""
Output consistency of simplify and VCF writing: strand agreement, self edges, determinism and inversion fields.
"""
import os
import subprocess
import sys
import warnings

import pytest

from pantree.graph import PangenomeGraph
from pantree.utils import node_complement, sequence_complement

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def load(fixture, ref_name, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return PangenomeGraph.from_gfa(os.path.join(DATA_DIR, fixture), ref_name=ref_name, **kwargs)


@pytest.mark.parametrize('fixture, ref_name, minimum_allele_length', [
    ('c4a_with_inversion_and_sequences.gfa', 'GRCh38', 50),
    ('simple_nested.gfa', 'ref', 1),
    ('self_edges_contraction.gfa', 'REF', 20),
])
def test_simplified_strands_are_reverse_complements(fixture, ref_name, minimum_allele_length):
    simplified = load(fixture, ref_name).simplify_subgraph(minimum_allele_length=minimum_allele_length)
    for node in simplified:
        assert simplified.nodes[node]['sequence'] == \
            sequence_complement(simplified.nodes[node_complement(node)]['sequence'])


def test_single_segment_is_written_forward(tmp_path):
    gfa = tmp_path / 'single.gfa'
    gfa.write_text('H\tVN:Z:1.1\nS\t1\tCAT\nW\tref\t0\tchr1\t0\t3\t>1\n')
    from pantree.simplify_graph import simplify_graph
    simplify_graph(str(gfa), str(tmp_path / 'out.gfa'), 'ref', 1, None, None)
    segments = [l.split('\t') for l in (tmp_path / 'out.gfa').read_text().splitlines() if l.startswith('S')]
    assert ['CAT'] == [s[2] for s in segments if s[3] == 'oi:J:["1"]']


def test_self_complementary_edge_is_a_variant(tmp_path):
    G = load('self_edges_genotype.gfa', 'REF')
    assert ('2_+', '2_-') in G.variant_edges
    G.write_vcf(os.path.join(DATA_DIR, 'self_edges_genotype.gfa'), str(tmp_path / 'out.vcf'), 'chr1')


def test_inversion_records_report_missing_ref_counts(tmp_path):
    G = load('inverted_reference_revisit.gfa', 'REF')
    vcf = tmp_path / 'out.vcf'
    G.write_vcf(os.path.join(DATA_DIR, 'inverted_reference_revisit.gfa'), str(vcf), 'chr1')
    records = [l.split('\t') for l in vcf.read_text().splitlines() if not l.startswith('#')]
    inversions = [r for r in records if 'VT=INV' in r[7]]
    assert inversions
    for r in inversions:
        assert 'RC=.;' in r[7] and 'AN=.;' in r[7]
        assert all(sample.split(':')[1] == '.' for sample in r[9:])


def test_simplify_output_does_not_depend_on_hash_seed(tmp_path):
    outputs = []
    for seed in ('1', '2'):
        out = tmp_path / f'out{seed}.gfa'
        script = ('from pantree.simplify_graph import simplify_graph; '
                  f'simplify_graph({os.path.join(DATA_DIR, "simple_nested.gfa")!r}, {str(out)!r}, "sample1", 50, '
                  'None, None)')
        subprocess.run([sys.executable, '-c', script], check=True, capture_output=True,
                       env={**os.environ, 'PYTHONHASHSEED': seed})
        outputs.append(out.read_text())
    assert outputs[0] == outputs[1]


def test_missing_reference_sample_warns():
    with pytest.warns(UserWarning, match='No walk of reference sample'):
        PangenomeGraph.from_gfa(os.path.join(DATA_DIR, 'simple_nested.gfa'), ref_name='not_a_sample')


def test_extra_reference_walks_are_reported_once():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        PangenomeGraph.from_gfa(os.path.join(DATA_DIR, 'c4a_with_inversion_and_sequences.gfa'), ref_name='HG00438')
    messages = [str(w.message) for w in caught if 'is used as the reference path' in str(w.message)]
    assert len(messages) == 1 and 'has 2 walks' in messages[0]


def test_sequence_complement_handles_soft_masking_and_ambiguity_codes():
    assert sequence_complement('ACGTacgtN') == 'NacgtACGT'
    assert sequence_complement('RYKMSWBDHV') == 'BDHVWSKMRY'
