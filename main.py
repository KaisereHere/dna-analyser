from visualize.plots import plot_gc_profile 
from tools.helpers import read_file
from analyse.sequence import read_fasta, open_read_frame

def main():
    data = read_file('data/rosalind_orf.txt')
    fasta = list(read_fasta(data).values())[0]
    fasta = 'ATGRER'
    protein = open_read_frame(fasta)
    print(protein)


if __name__ == "__main__":
    main()