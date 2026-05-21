import pytest

from analyse.sequence import reversed_complement, transcribe, count_nucleotides, read_fasta, calculate_gc_content, calculate_hamming_distance
from tools.helpers import read_file

class TestSequence:

    strand_test = 'AAAACCCGGT'
    fasta_test = '>fd3\nCT\nGA\n>fs3\nGTC'

    def test_reversed_completement(self):
        assert 'ACCGGGTTTT' == reversed_complement(self.strand_test)

    def test_reversed_completement_empty(self):
        assert '' == reversed_complement('')

    def test_transcribe_empty(self):
        assert transcribe('') == ''

    def test_transcribe_lowercase(self):
        assert transcribe('aaaacccggt') == 'AAAACCCGGU'

    def test_transcribe(self):
        assert transcribe(self.strand_test) == 'AAAACCCGGU'

    def test_count_nucleotides_empty(self):  
        assert count_nucleotides('') == {'A':0, 'C':0, 'G':0, 'T':0}

    def test_count_nucleotides_random(self):  
        assert count_nucleotides('dddppbbbb343433434f') == {'A':0, 'C':0, 'G':0, 'T':0}

    def test_count_nucleotides(self):  
        assert count_nucleotides(self.strand_test) == {'A':4, 'C':3, 'G':2, 'T':1}

    def test_read_file(self):
        test_bcode = '>Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC\nTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT\nATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC\nTGGGAACCTGCGGGCAGTAGGTGGAAT'
        data = read_file('data/test_fasta.fasta')
        assert data == test_bcode

    def test_read_fasta(self):
        assert read_fasta(self.fasta_test) == {'fd3':'CTGA', 'fs3':'GTC'}
    
    def test_read_fasta_empty(self):
        assert read_fasta('') == {}

    def test_calculate_gc_content_empty(self):
        assert calculate_gc_content('') == 0

    def test_calculate_gc_content(self):   
        strand = 'CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT'
        assert abs(calculate_gc_content(strand) - 60.919540) <= 0.001

    def test_calculate_hamming_distance(self):
        test_strand_1 = 'GAGCCTACTAACGGGAT'
        test_strand_2 = 'CATCGTAATGACGGCCT'
        assert calculate_hamming_distance(test_strand_1, test_strand_2) == 7

    def test_calculate_hamming_distance_empty(self):
        test_strand_1 = ''
        test_strand_2 = ''
        assert calculate_hamming_distance(test_strand_1, test_strand_2) == 0

    def test_calculate_hamming_distance_length(self):
        with pytest.raises(ValueError):
            calculate_hamming_distance('A', 'AT')