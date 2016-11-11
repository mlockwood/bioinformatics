__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from class_II.string_composition import *


def test_string_composition():
    assert string_composition(5, 'CAATCCAAC') == 'AATCC\nATCCA\nCAATC\nCCAAC\nTCCAA'


def test_string_reconstruction_pre_ordered():
    assert string_reconstruction_pre_ordered('ACCGA\nCCGAA\nCGAAG\nGAAGC\nAAGCT') == 'ACCGAAGCT'


def test_get_prefix():
    assert get_prefix('TAA') == 'TA'


def test_get_suffix():
    assert get_suffix('TAA') == 'AA'


def test_get_overlap_graph():
    assert get_overlap_graph('ATGCG\nGCATG\nCATGC\nAGGCA\nGGCAT'
                             ) == 'AGGCA -> GGCAT\nCATGC -> ATGCG\nGCATG -> CATGC\nGGCAT -> GCATG\n'


def test_de_bruijn_graph_from_string():
    assert de_bruijn_graph_from_string(4, 'AAGATTCTCTAAGA') == (
        'AAG -> AGA,AGA\nAGA -> GAT\nATT -> TTC\nCTA -> TAA\nCTC -> TCT\nGAT -> ATT\nTAA -> AAG\nTCT -> CTA,CTC\n' +
        'TTC -> TCT\n'
    )
