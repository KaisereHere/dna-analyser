from visualize.plots import plot_gc_profile 
from tools.helpers import read_file
from analyse.sequence import read_fasta, translate, transcribe, rna_splicing

def main():
    raw_strand_data = read_file('data/rosalind_spl_or.txt')
    data = read_file('data/rosalind_splc.txt')
    raw_strand = list(read_fasta(raw_strand_data).values())[0]
    introns = read_fasta(data).values()
    strand = rna_splicing(raw_strand, introns)

    print(translate(transcribe(strand)))

if __name__ == "__main__":
    main()