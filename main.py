from analyse.sequence import transcribe, count_nucleotides, reversed_complement, read_fasta, calculate_gc_content
from tools.helpers import read_file

def main():
    fasta = read_fasta(read_file('data/rosalind_gc.txt'))
    biggest_gc = 0
    name_strand = ''
    for name, strand in fasta.items():
        if calculate_gc_content(strand) > biggest_gc:
            biggest_gc = calculate_gc_content(strand)
            name_strand = name
    print(name_strand, biggest_gc)

if __name__ == "__main__":
    main()