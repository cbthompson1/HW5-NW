# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # Initialize alignment and gap matrices
        self._align_matrix = np.zeros([len(seqA) + 1, len(seqB) + 1])
        self._gapA_matrix = np.zeros([len(seqA) + 1, len(seqB) + 1])
        self._gapB_matrix = np.zeros([len(seqA) + 1, len(seqB) + 1])
        for i in range(self._align_matrix.shape[0]):
            for j in range(self._align_matrix.shape[1]):
                if i == 0:
                    # Case: SeqA edge initialization to -inf
                    self._gapA_matrix[i][j] = -np.inf
                    self._align_matrix[i][j] = -np.inf
                else:
                    # Case: Calculate gap creation/extension score along seqA.
                    self._gapA_matrix[i][j] = max(
                        self._align_matrix[i-1][j] + self.gap_open + self.gap_extend,
                        self._gapA_matrix[i-1][j] + self.gap_extend
                    )
                if j == 0:
                    # Case: SeqB edge initialization to -inf
                    self._gapB_matrix[i][j] = -np.inf
                    self._align_matrix[i][j] = -np.inf
                else:
                    # Case: Calculate gap creation/extension score along seqB.
                    self._gapB_matrix[i][j] = max(
                        self._align_matrix[i][j-1] + self.gap_open + self.gap_extend,
                        self._gapB_matrix[i][j-1] + self.gap_extend
                    )
                if i == 0 and j == 0:
                    # Special Case: (0,0) should be 0 as an end state.
                    self._align_matrix[i][j] = 0
                elif i >= 1 and j >= 1:
                    # Case: Diagonal calculation of align matrix.
                    # Use either the sub_dict or default to +1/-1.
                    key = (seqA[i-1], seqB[j-1])
                    if key in self.sub_dict:
                        score = self.sub_dict[key]
                    else:
                        print(f"WARNING: key pairing {key} not in sub dict")
                        score = 1 if seqA[i-1] == seqB[j-1] else -1
                    self._align_matrix[i][j] = max(
                        self._align_matrix[i-1][j-1] + score,
                        self._gapA_matrix[i-1][j-1] + score,
                        self._gapB_matrix[i-1][j-1] + score,
                    )
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """

        # Easy first step, just find the max align value to report at the end.
        self.alignment_score = max(
            self._align_matrix[-1][-1],
            self._gapA_matrix[-1][-1],
            self._gapB_matrix[-1][-1]
        )

        # Start at the bottom right and continue until the top left is reached.
        matrix_index = (len(self._seqA), len(self._seqB))
        while matrix_index != (0,0):
            # Sequence values are one less than the matrix index because of the
            # additional matrix row/col.
            a_index = matrix_index[0] - 1
            b_index = matrix_index[1] - 1
            
            # Look at the three values in the current cell.
            align_val = self._align_matrix[matrix_index[0]][matrix_index[1]]
            gap_a_val = self._gapA_matrix[matrix_index[0]][matrix_index[1]]
            gap_b_val = self._gapB_matrix[matrix_index[0]][matrix_index[1]]

            # Take the maximum value as indication of what direction to move.
            # Also, as the function progresses, append to the sequences.
            if align_val >= max(gap_a_val, gap_b_val):
                # Case align: move diagonal, match/mismatch.
                self.seqA_align = self._seqA[a_index] + self.seqA_align
                self.seqB_align = self._seqB[b_index] + self.seqB_align
                matrix_index = (matrix_index[0] - 1, matrix_index[1] - 1)
            elif gap_a_val >= max(gap_b_val, align_val):
                # Case gap (A): make a gap for A sequence, move up a row.
                self.seqA_align = self._seqA[a_index] + self.seqA_align
                self.seqB_align = '-' + self.seqB_align
                matrix_index = (matrix_index[0] - 1, matrix_index[1])
            elif gap_b_val >= max(gap_a_val, align_val):
                # Case gap (B): make a gap for B sequence, move down a col.
                self.seqA_align = '-' + self.seqA_align
                self.seqB_align = self._seqB[b_index] + self.seqB_align
                matrix_index = (matrix_index[0], matrix_index[1] - 1)
        # Return the backtracing work.
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
