from analyse.models import UniProtRecord, Sequence, DNASequence, ProteinSequence, FastaSequence
from tools.reports import get_info_protein
from analyse.sequence import prosite_parse, read_fasta, transition_to_transversion_ratio, align
from tools.helpers import uni_prot_metadata_parse, get_protein_uniprot, read_file

def main():
 #   get_info_protein('P68871')
#    get_info_protein('P69905')

    hbb = UniProtRecord(get_protein_uniprot('P68871'))
    hba = UniProtRecord(get_protein_uniprot('P69905'))

    data = read_file('data/rosalind_cons.txt')
    print(FastaSequence.array(data))

if __name__ == "__main__":
    main()