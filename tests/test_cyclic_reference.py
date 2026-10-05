"""
Reference walks that revisit a node (tandem or inverted duplications on the reference itself).
"""
import os
import warnings

import pytest

from pantree.dfs import dfs_methods
from pantree.graph import PangenomeGraph

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def assert_single_rooted_spanning_tree(G):
    tree = G.reference_tree
    roots = [node for node in tree if tree.in_degree(node) == 0]
    assert roots == ['+_terminus_+']
    assert max(degree for _, degree in tree.in_degree()) == 1
    positive_nodes = [node for node, direction in G.nodes(data='direction') if direction == 1]
    assert all(node in tree for node in positive_nodes)


def test_reference_revisits_node_in_opposite_orientation():
    gfa = os.path.join(DATA_DIR, 'inverted_reference_revisit.gfa')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        G = PangenomeGraph.from_gfa(gfa, ref_name='REF')
    assert_single_rooted_spanning_tree(G)
    assert G.variant_edges == {('3_+', '2_-'), ('2_-', '4_+')}
    assert all(G.is_inversion(edge) for edge in G.variant_edges)


@pytest.mark.parametrize('ref_name', ['HG02257', 'HG03540'])
def test_reference_revisits_node_in_same_orientation(ref_name):
    gfa = os.path.join(DATA_DIR, 'c4a_with_inversion_and_sequences.gfa')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        G = PangenomeGraph.from_gfa(gfa, ref_name=ref_name)
    assert_single_rooted_spanning_tree(G)


def test_simplify_keeps_cyclic_reference_walk_connected():
    # Segment 3 is only used by the reference, so the hairpin's inversion edges are its only links onward
    gfa = os.path.join(DATA_DIR, 'reference_only_hairpin.gfa')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        G = PangenomeGraph.from_gfa(gfa, ref_name='REF')
    pruned = G.delete_small_variants(minimum_allele_length=50)
    for edge in zip(G.reference_path[:-1], G.reference_path[1:]):
        if G.has_edge(*edge):
            assert pruned.has_edge(*edge)
    G.simplify_subgraph(minimum_allele_length=50)


@pytest.mark.parametrize('fixture, minimum_allele_length', [
    # cycle made only of tree and cross edges unless the tree is a genuine DFS tree (topological sort)
    ('cyclic_reference_cross_edges.gfa', 11),
    # part of the graph entered only through inversion edges, unreachable from the + terminus
    ('cyclic_reference_unreachable.gfa', 29),
    # path contracted together with its own complement in simplify
    ('cyclic_reference_hairpin_contraction.gfa', 30),
    ('inverted_reference_revisit.gfa', 1),
    ('reference_only_hairpin.gfa', 50),
])
@pytest.mark.parametrize('dfs_method', ['max_weight', 'contiguous'])
def test_cyclic_reference_end_to_end(fixture, minimum_allele_length, dfs_method, tmp_path):
    gfa = os.path.join(DATA_DIR, fixture)
    priority_dict = {'REF': 0, 'ALT': 1} if dfs_method == 'contiguous' else None
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        G = PangenomeGraph.from_gfa(gfa, ref_name='REF', dfs_method=dfs_methods[dfs_method],
                                    priority_dict=priority_dict)
    assert_single_rooted_spanning_tree(G)
    G.simplify_subgraph(minimum_allele_length=minimum_allele_length)
    G.write_vcf(gfa, str(tmp_path / 'out.vcf'), 'chr1')


def test_acyclic_reference_does_not_warn():
    gfa = os.path.join(DATA_DIR, 'c4a_with_inversion_and_sequences.gfa')
    with warnings.catch_warnings():
        warnings.filterwarnings('error', message='Reference walk is cyclic')
        G = PangenomeGraph.from_gfa(gfa, ref_name='GRCh38')
    assert_single_rooted_spanning_tree(G)


def test_acyclic_reference_with_unreachable_region(tmp_path):
    # ALT enters segment 4 only through self edges, so nothing reaches 4_+ from the + terminus
    gfa = os.path.join(DATA_DIR, 'acyclic_reference_unreachable.gfa')
    with pytest.warns(UserWarning, match='unreachable'):
        G = PangenomeGraph.from_gfa(gfa, ref_name='REF')
    assert_single_rooted_spanning_tree(G)
    G.write_vcf(gfa, str(tmp_path / 'out.vcf'), 'chr1')


def test_revisited_reference_nodes_keep_first_position():
    # REF walk >1>2>3<2>4 with 8 bp segments: the revisit of segment 2 does not move it
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        G = PangenomeGraph.from_gfa(os.path.join(DATA_DIR, 'inverted_reference_revisit.gfa'), ref_name='REF')
    assert G.position('2_+') == 16
    assert G.position('4_+') == 40
