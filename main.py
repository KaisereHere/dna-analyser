from analyse.sequence import transcribe, count_nucleotides, reversed_complement, read_fasta, calculate_gc_content, calculate_hamming_distance
from tools.helpers import read_file

def main():
    strand1, strand2 = read_file('data/rosalind_hamm.txt').split('\n')[:2]
    print(calculate_hamming_distance(strand1, strand2))


if __name__ == "__main__":
    main()