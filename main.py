from analyse.sequence import read_fasta, spliced_motif, shared_spliced_motif, edit_distance
from tools.helpers import read_file


def main():
    data = read_file('data/rosalind_edit.txt')

    fasta = list(read_fasta(data).values())

    res = edit_distance(fasta[0], fasta[1])

    print(res)
    #fasta = list(read_fasta(data).values())
    #print(shared_spliced_motif(fasta[0], fasta[1]))
if __name__ == "__main__":
    main()
