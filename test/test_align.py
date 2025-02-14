# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Unit test for Needleman Wunsch align function.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    need_wun = NeedlemanWunsch(
        'substitution_matrices/BLOSUM62.mat',
        gap_open=-10,
        gap_extend=-1
    )
    score, seq1_align, seq2_align = need_wun.align(seq1, seq2)
    
    # 3 matching (+3), one gap creation (-10) of size 1 (-1).
    assert score == 4

    # Seq1 should remain the same, seq2 should have one dash.
    assert seq1 == seq1_align
    assert len(seq2_align) == 4
    assert '-' in seq2_align
    
    # These are what the calculated matrices should look like.
    assert np.allclose(need_wun._align_matrix, np.array([
        [  0., -np.inf, -np.inf, -np.inf],
        [-np.inf,   5., -11., -13.],
        [-np.inf, -12.,   4.,  -8.],
        [-np.inf, -12.,  -1.,   5.],
        [-np.inf, -14.,  -6.,   4.]
    ]))
    assert np.allclose(need_wun._gapA_matrix, np.array([
        [-np.inf, -np.inf, -np.inf, -np.inf],
        [-11., -np.inf, -np.inf, -np.inf],
        [-12.,  -6., -22., -24.],
        [-13.,  -7.,  -7., -19.],
        [-14.,  -8.,  -8.,  -6.]
    ]))
    assert np.allclose(need_wun._gapB_matrix, np.array([
        [-np.inf, -11., -12., -13.],
        [-np.inf, -np.inf,  -6.,  -7.],
        [-np.inf, -np.inf, -23.,  -7.],
        [-np.inf, -np.inf, -23., -12.],
        [-np.inf, -np.inf, -25., -17.]
    ]))

def test_nw_backtrace():
    """
    Confirm backtrace works.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    need_wun = NeedlemanWunsch(
        'substitution_matrices/BLOSUM62.mat',
        gap_open=-10,
        gap_extend=-1
    )
    score, seqA, seqB = need_wun.align(
        seq3,
        seq4
    )

    # From README, should be 17 and have the following aligned sequences.
    assert score == 17
    assert seqA == 'MAVHQLIRRP'
    assert seqB == 'M---QLIRHP'