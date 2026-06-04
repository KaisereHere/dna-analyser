import pytest

from analyse.sequence import reversed_complement, transcribe, count_nucleotides, read_fasta, calculate_gc_content, calculate_hamming_distance
from analyse.sequence import find_motif, translate, gc_sliding_window, dominant_probability, find_monoisotopic_mass, strand_profile, consensus_strand
from analyse.sequence import open_read_frame
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

    def test_find_motif(self):
        assert find_motif('GATATATGCATATACTT', 'ATAT') == [2, 4, 10]

    def test_find_motif_empty(self):
        assert find_motif('', '') == []

    def test_find_motif_overload(self):
        assert find_motif('A', 'AA') == []

    def test_translate(self):
        assert translate('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA') == 'MAMAPRTEINSTRING'

    def test_translate_empty(self):
        assert translate('') == ''

    def test_translate_false_codon(self):
        with pytest.raises(ValueError):
            translate('UUT')

    def test_translate_false_nucleotide_count(self):
        with pytest.raises(ValueError):
            translate('U')

    def test_gc_sliding_window(self):
        
        gcs_list = gc_sliding_window('CCGAACGGA', 3) 

        for gc, gc_test in zip(gcs_list, [100, 66.666, 33.333, 33.333, 66.666, 100.0, 66.666]):
            assert int(gc) == int(gc_test)
        
    def test_gc_window_too_large(self):
        with pytest.raises(IndexError):
            gc_sliding_window('CCG', 5)

    def test_gc_window_zero(self):
        with pytest.raises(IndexError):
            gc_sliding_window('CCG', 0)


    def test_dominant_probability(self):
        assert abs(dominant_probability(2,2,2) - 0.78333) <= 0.00001

    def test_dominant_probability_negative(self):
        with pytest.raises(ValueError):
            dominant_probability(-1, -1, -1)

    
    def test_dominant_probability_int(self):
        with pytest.raises(TypeError):
            dominant_probability(1.5, 0.5, 2)

    def test_dominant_probability_zero(self):
        assert dominant_probability(0, 0, 0) == 0

    def test_find_monoisotopic_mass(self):
        assert abs(find_monoisotopic_mass('SKADYEK')-821.392) <= 0.0001

    def test_find_monoisotopic_mass_empty(self):
        assert find_monoisotopic_mass('') == 0

    def test_find_monoisotopic_mass_random_symbols(self):
        with pytest.raises(ValueError):
            find_monoisotopic_mass('UB')

    def test_strand_profile(self):
        test_set = {'Rosalind_1': 'ATCCAGCT', 'Rosalind_2': 'GGGCAACT', 'Rosalind_3': 'ATGGATCT', 'Rosalind_4': 'AAGCAACC', 'Rosalind_5': 'TTGGAACT', 'Rosalind_6': 'ATGCCATT', 'Rosalind_7': 'ATGGCACT'}
        test_profile = {
            'A': [5, 1, 0, 0, 5, 5, 0, 0],
            'C': [0, 0, 1, 4, 2, 0, 6, 1],
            'G': [1, 1, 6, 3, 0, 1, 0, 0],
            'T': [1, 5, 0, 0, 0, 1, 1, 6]
        }
        assert strand_profile(test_set) == test_profile

    def test_strand_profile_wrong_size(self):
        with pytest.raises(ValueError):
            strand_profile({'a': "AC", 'b': "GGG"})


    def test_strand_profile_wrong_nucleotide(self):
        with pytest.raises(ValueError):
            strand_profile({"a":"R"})

    def test_strand_profile_empty(self):
        with pytest.raises(ValueError):
            strand_profile({})
    
    def test_consensus_strand(self):
        assert consensus_strand({'A': [5, 1, 0, 0, 5, 5, 0, 0], 'G': [1, 1, 6, 3, 0, 1, 0, 0], 'C': [0, 0, 1, 4, 2, 0, 6, 1], 'T': [1, 5, 0, 0, 0, 1, 1, 6]}) == 'ATGCAACT'

    def test_consensus_strand_empty(self):
        assert consensus_strand({'A':[],'G':[], 'C':[],'T':[]}) == '' 

    def test_open_read_frame(self):
        test_strand = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
        assert set(['MLLGSFRLIPKETLIQVAGSSPCNLS', 'M', 'MGMTPRLGLESLLE', 'MTPRLGLESLLE']) == open_read_frame(test_strand)

    def test_open_read_frame_empty(self):
        assert open_read_frame('') == set()
    
    def test_open_read_random(self):
        with pytest.raises(KeyError):
            open_read_frame('ATGRE')