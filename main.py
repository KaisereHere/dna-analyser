from analyse.models import UniProtRecord, Sequence, DNASequence, ProteinSequence, FastaSequence
from tools.reports import protein_report, compare_2_proteins_report

from visualize.plots import sliding_window_plot
from analyse.sequence import prosite_parse, read_fasta, transition_to_transversion_ratio, align, identity_sliding_window
from tools.helpers import uni_prot_metadata_parse, get_protein_uniprot, read_file

def main():
    data = read_file('data/rosalind_long.txt')
    fasta = FastaSequence.array(data)
    print(DNASequence.assemble(fasta))
    
if __name__ == "__main__":
    main()