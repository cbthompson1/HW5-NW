# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import re

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """

    order = []
    filepaths = [
        "./data/Gallus_gallus_BRD2.fa",
        "./data/Mus_musculus_BRD2.fa",
        "./data/Balaeniceps_rex_BRD2.fa",
        "./data/tursiops_truncatus_BRD2.fa"
    ]
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    needleman_wunsch = NeedlemanWunsch(
        'substitution_matrices/BLOSUM62.mat',
        gap_open=-10,
        gap_extend=-1
    )
    for fp in filepaths:
        comparison_seq, header = read_fasta(fp)
        score, hs_align, comp_align = needleman_wunsch.align(hs_seq, comparison_seq)
        order.append((score, header, hs_align, comp_align))
    order.sort(reverse=True)

    print("Species in order of alignment with Homo Sapiens and their scores:")
    for seq in order:
        score, header, hs_align, comp_align = seq
        comp_name = _get_species_name(header)
        print(f'{comp_name}: {score}')

def _get_species_name(header):
    """Regex search the species name out."""
    return re.search(r'OS=(.+) OX=', header).group(1)

if __name__ == "__main__":
    main()
