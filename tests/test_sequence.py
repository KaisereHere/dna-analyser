from analyse.sequence import reversed_complement, transcribe, count_nucleotides, read_fasta, calculate_gc_content
from tools.helpers import read_file

class TestSequence:

    strand_test = 'AAAACCCGGT'

    test_fasta = {
            "Rosalind_6404": 'CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC' +
            'TCCCACTAATAATTCTGAGG',
            "Rosalind_5959": 'CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT' +
            'ATATCCATTTGTCAGCAGACACGC',
            "Rosalind_0808": 'CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC' +
            'TGGGAACCTGCGGGCAGTAGGTGGAAT'
    }

    def test_reversed_completement(self):
        assert 'ACCGGGTTTT' == reversed_complement(self.strand_test)

    def test_transcribe(self):
        assert transcribe(self.strand_test) == 'AAAACCCGGU'

    def test_count_nucleotides(self):  
        assert count_nucleotides(self.strand_test) == {'A':4, 'C':3, 'G':2, 'T':1}

    def test_read_fasta(self):
        data = read_file('data/test_fasta.fasta')
        assert read_fasta(data) == self.test_fasta

    def test_calculate_gc_content(self):   
        assert abs(calculate_gc_content(self.test_fasta['Rosalind_0808']) - 60.919540) <= 0.001
