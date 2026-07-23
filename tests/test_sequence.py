import pytest

from analyse.sequence import reversed_complement, transcribe, count_nucleotides, read_fasta, calculate_gc_content, calculate_hamming_distance
from analyse.sequence import find_motif, translate, gc_sliding_window, dominant_probability, find_monoisotopic_mass, strand_profile, consensus_strand
from analyse.sequence import open_read_frame, rna_splicing, dominant_offspring, find_shared_motif, independent_alleles, find_component, find_components
from analyse.sequence import k_mer_sliding_window, assembly_DNA, inferring_mRNA_possablilities, de_bruijn, _get_motif_profile
from analyse.sequence import spliced_motif, shared_spliced_motif, edit_distance, overlapping_sequences, enumerate_k_mers, build_adjacency_list
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
        with pytest.raises(ValueError):
            gc_sliding_window('CCG', 5)

    def test_gc_window_zero(self):
        with pytest.raises(ValueError):
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
    

    def test_rna_splicing(self):
        assert rna_splicing('ATGGCAC', ['GG', 'AT']) == 'CAC'

    def test_rna_splicing_too_large(self):
        with pytest.raises(ValueError):
            rna_splicing('A', ['AA'])

    def test_rna_splicing_empty(self):
        assert rna_splicing('', []) == ''

    def test_dominant_offspring(self):
        assert dominant_offspring(1, 0, 0, 1, 0, 1) == 3.5

    def test_find_shared_motif(self):
        test_data = {'Rosalind_1': 'GATTACA', 'Rosalind_2': 'TAGACCA', 'Rosalind_3': 'ATACA'}
        result = find_shared_motif(test_data)
        assert result in ("AC", "CA", "TA")

    def test_find_shared_motif_empty(self):
        with pytest.raises(ValueError):
            find_shared_motif({})

    def test_find_shared_motif_no_shared(self):
        assert find_shared_motif({'d': 'AC', 'b': "GT"}) == ''

    def test_independent_alleles(self):
        assert abs(independent_alleles(2, 1)-0.684) <= 0.001

    def test_independent_alleles_zeros(self):
        assert independent_alleles(0, 0) == 1

    def test_independent_alleles_more_expected(self):
        assert independent_alleles(1, 5) == 0 

    def test_spliced_motif(self):
        sequence = 'ACGTACGTGACG'
        subsequence = 'GTA' 
        indices = spliced_motif(sequence, subsequence)

        for index, value in enumerate(indices):
            assert subsequence[index] == sequence[value] 

    def test_spliced_motif_empty_subsequence(self):
        with pytest.raises(ValueError):
            spliced_motif('GGA', '')

    def test_spliced_motif_not_found(self):
        assert spliced_motif("AGG", "T") == []

    def test_find_shared_motif_border(self):
        assert find_shared_motif({"a":"AGGG","b": "CCCA"}) == "A"

    def test_shared_spliced_motif(self):
        assert ''.join(shared_spliced_motif('ACACTGTGA', 'AACCTTGG')) in('AACTGG', 'ACCTGG') 

    def test_shared_spliced_motif_similar(self):
        assert ''.join(shared_spliced_motif('AA', 'AA')) == 'AA'   

    def test_shared_spliced_motif_same_letters(self):
        assert shared_spliced_motif("AAAA", "AA") == 'AA'

    def test_shared_spliced_motif_no_match(self):
        assert shared_spliced_motif('AAA', 'BB') == ''

    def test_shared_spliced_motif_empty(self):
        assert shared_spliced_motif('', 'AAA') == ''

    def test_edit_distance(self):
        assert edit_distance('PLEASANTLY', "MEANLY") == 5
        
    def test_edit_distance_empty(self):
        assert edit_distance('', '') == 0

    def test_edit_distance_same_strings(self):
        assert edit_distance('AA', 'AA') == 0 

    def test_edit_distance_second_larger(self):
        assert edit_distance('A', 'B') == 1

    def test_edit_distance_three(self):
        assert edit_distance('', 'ABC') == 3

    def test_overlapping_sequences(self):
        sequences = {
            'Rosalind_0498': 'AAATAAA',
            'Rosalind_2391': 'AAATTTT',
            'Rosalind_2323': 'TTTTCCC',
            'Rosalind_0442': 'AAATCCC',
            'Rosalind_5013': 'GGGTGGG'
        }

        ordered = [('Rosalind_0498', 'Rosalind_2391'),
                   ('Rosalind_0498', 'Rosalind_0442'),
                   ('Rosalind_2391', 'Rosalind_2323')
        ]

       
        assert overlapping_sequences(sequences) == set(ordered)
    

    def test_overlapping_sequences_empty(self):
        assert overlapping_sequences({}) == set()

     
    def test_overlapping_sequences_k_large(self):
        with pytest.raises(ValueError):
            overlapping_sequences({'a': 'AAA'}, 5) 

    def test_enumerate_k_mers(self):
        answer = [
            'AA',
            'AC',
            'AG',
            'AT',
            'CA',
            'CC',
            'CG',
            'CT',
            'GA',
            'GC',
            'GG',
            'GT',
            'TA',
            'TC',
            'TG',
            'TT'
]
    
        assert enumerate_k_mers(['A', 'C', 'G', 'T'], 2) == answer

    def test_enumerate_k_mers_empty(self):
        assert enumerate_k_mers(['A', 'C'], 0) == ['']

    def test_build_adjacency_list(self):
        assert build_adjacency_list([1,2], [(1, 2)]) == {1: set([2]), 2: set([1])}

    def test_build_adjacency_list_isolated_vertex(self):
        assert build_adjacency_list([1], []) == {1: set()}

    def test_find_component(self):
        assert set(find_component(2, {2: set([5, 1]), 1: set([2, 5]), 5: set([2, 1]), 10: set([11]), 11:set([10])})) == set([1, 5, 2 ])

    def test_find_component_empty(self):
        assert find_component(1, {}) == []

    def test_find_components(self):
        assert find_components([1, 2], [(1, 2)]) in ([[1, 2]], [[2, 1]])

    def test_k_mer_sliding_window(self):
        seq = 'CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGGCCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGTTTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCAAATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCGGGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGACTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTACCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG'
        res ={'AAAA': 4, 'AAAC': 1, 'AAAG': 4, 'AAAT': 3, 'AACA': 0, 'AACC': 1, 'AACG': 1, 'AACT': 5, 'AAGA': 1, 'AAGC': 3, 'AAGG': 1, 'AAGT': 2, 'AATA': 2, 'AATC': 1, 'AATG': 2, 'AATT': 0, 'ACAA': 1, 'ACAC': 1, 'ACAG':3, 'ACAT': 1, 'ACCA': 2, 'ACCC': 1, 'ACCG': 3, 'ACCT': 1, 'ACGA': 1, 'ACGC': 1, 'ACGG': 1, 'ACGT': 2, 'ACTA': 2, 'ACTC': 5, 'ACTG': 1, 'ACTT': 3, 'AGAA': 0, 'AGAC': 2, 'AGAG': 2, 'AGAT': 1, 'AGCA': 1, 'AGCC': 1, 'AGCG': 1, 'AGCT': 3, 'AGGA': 1, 'AGGC': 0, 'AGGG': 0, 'AGGT': 1, 'AGTA': 5, 'AGTC': 5, 'AGTG': 1, 'AGTT': 5, 'ATAA': 0, 'ATAC': 2, 'ATAG': 0, 'ATAT': 2, 'ATCA': 1, 'ATCC': 2, 'ATCG': 1, 'ATCT': 1, 'ATGA': 1, 'ATGC': 2, 'ATGG': 0, 'ATGT': 1, 'ATTA': 0, 'ATTC': 0, 'ATTG': 1, 'ATTT': 1, 'CAAA': 3, 'CAAC': 2, 'CAAG': 1, 'CAAT': 0, 'CACA': 3, 'CACC': 2, 'CACG': 3, 'CACT': 0, 'CAGA': 0, 'CAGC': 2, 'CAGG': 0, 'CAGT': 8, 'CATA': 0, 'CATC': 0, 'CATG': 1, 'CATT': 0, 'CCAA': 2, 'CCAC': 1, 'CCAG': 3, 'CCAT': 0, 'CCCA': 0, 'CCCC': 0, 'CCCG': 1, 'CCCT': 4, 'CCGA': 3, 'CCGC': 2, 'CCGG': 1, 'CCGT': 1, 'CCTA': 3, 'CCTC': 1,'CCTG': 2, 'CCTT': 1, 'CGAA': 3, 'CGAC': 1, 'CGAG': 2, 'CGAT': 1, 'CGCA': 2, 'CGCC': 1, 'CGCG': 1, 'CGCT': 1, 'CGGA': 2, 'CGGC': 3, 'CGGG': 2, 'CGGT': 1, 'CGTA': 1, 'CGTC': 0, 'CGTG': 1, 'CGTT': 1, 'CTAA': 3, 'CTAC': 2, 'CTAG': 1, 'CTAT': 2, 'CTCA': 6, 'CTCC': 2, 'CTCG': 1, 'CTCT': 1, 'CTGA': 1, 'CTGC': 2, 'CTGG': 3, 'CTGT': 3, 'CTTA': 3, 'CTTC': 2, 'CTTG': 3, 'CTTT': 0, 'GAAA': 3, 'GAAC': 2, 'GAAG': 1, 'GAAT': 1, 'GACA': 0, 'GACC': 0, 'GACG': 1, 'GACT': 4, 'GAGA': 3, 'GAGC': 0, 'GAGG': 1, 'GAGT': 5, 'GATA': 0, 'GATC': 2, 'GATG': 0, 'GATT': 1, 'GCAA': 2, 'GCAC': 1, 'GCAG': 3, 'GCAT': 0, 'GCCA': 1, 'GCCC': 2, 'GCCG': 2, 'GCCT': 1, 'GCGA': 1, 'GCGC': 0, 'GCGG': 3, 'GCGT': 0, 'GCTA': 0, 'GCTC': 4, 'GCTG': 5, 'GCTT': 0, 'GGAA': 3, 'GGAC': 0, 'GGAG': 2, 'GGAT': 1, 'GGCA': 1, 'GGCC': 3, 'GGCG': 0, 'GGCT': 3, 'GGGA': 2, 'GGGC': 2, 'GGGG': 1, 'GGGT': 1, 'GGTA': 0, 'GGTC': 2, 'GGTG': 1, 'GGTT': 0, 'GTAA': 2, 'GTAC': 2, 'GTAG': 1, 'GTAT': 2, 'GTCA': 0, 'GTCC': 2, 'GTCG': 2, 'GTCT': 5, 'GTGA': 2, 'GTGC': 2, 'GTGG': 1, 'GTGT': 1, 'GTTA': 2, 'GTTC': 1, 'GTTG': 2, 'GTTT': 2, 'TAAA': 2, 'TAAC': 2, 'TAAG': 1, 'TAAT': 1, 'TACA': 3, 'TACC': 4, 'TACG': 0, 'TACT': 2, 'TAGA': 1, 'TAGC': 1, 'TAGG': 0, 'TAGT': 1, 'TATA': 2, 'TATC': 2, 'TATG': 1, 'TATT': 1, 'TCAA': 1, 'TCAC': 5, 'TCAG': 2, 'TCAT': 0, 'TCCA': 3, 'TCCC': 2, 'TCCG': 1, 'TCCT': 1, 'TCGA': 2, 'TCGC': 2, 'TCGG': 3, 'TCGT': 0, 'TCTA': 3, 'TCTC': 0, 'TCTG': 1, 'TCTT': 3, 'TGAA': 1, 'TGAC':2, 'TGAG': 3, 'TGAT': 0, 'TGCA': 2, 'TGCC': 1, 'TGCG': 2, 'TGCT': 2, 'TGGA': 1, 'TGGC': 2, 'TGGG': 3, 'TGGT': 0, 'TGTA': 1, 'TGTC': 2, 'TGTG': 3, 'TGTT': 1, 'TTAA': 1, 'TTAC': 3, 'TTAG': 1, 'TTAT': 0, 'TTCA': 1, 'TTCC': 1, 'TTCG': 3, 'TTCT': 0, 'TTGA': 2, 'TTGC': 1, 'TTGG': 2, 'TTGT': 2, 'TTTA': 0, 'TTTC': 2, 'TTTG': 1, 'TTTT': 1} 
        assert k_mer_sliding_window(seq, 4) == res


    def test_k_mer_sliding_window_index_error(self):
        with pytest.raises(ValueError): 
            k_mer_sliding_window('A',2)

    
    def test_k_mer_sliding_window_zero(self):
        with pytest.raises(ValueError):
            k_mer_sliding_window('CCG', 0)

    def test_assembly_DNA(self):
        assert assembly_DNA({'Rosalind_56': 'ATTAGACCTG', 'Rosalind_57': 'CCTGCCGGAA', 'Rosalind_58': 'AGACCTGCCG', 'Rosalind_59': 'GCCGGAATAC'}) == 'ATTAGACCTGCCGGAATAC'
    
    def test_inferring_mRNA_possabilities(self):
        assert inferring_mRNA_possablilities('MA') == 12

    def test_inferring_mRNA_possabilities_zero(self):
        assert inferring_mRNA_possablilities('') == 3

    def test_de_bruijn(self):
        res = set([
            ('ATC', 'TCA'),
            ('ATG', 'TGA'),
            ('ATG', 'TGC'), 
            ('CAT', 'ATC'),
            ('CAT', 'ATG'),
            ('GAT', 'ATG'),
            ('GCA', 'CAT'),
            ('TCA', 'CAT'),
            ('TGA', 'GAT')
])
        
        inp = [
        'TGAT',
        'CATG',
        'TCAT',
        'ATGC',
        'CATC',
        'CATC'
]
        assert de_bruijn(inp) == res

    def test_de_bruijn_empty(self):
        assert de_bruijn([]) == set()

    def test_de_bruijn_single_symbol(self):
        assert de_bruijn(['A']) == set([('', '')])

    def test__get_motif_profile(self):  
        res = {0: ['N'], 1: ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'S', 'T', 'W', 'Y', 'V'], 2: ['S', 'T'], 3: ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'S', 'T', 'W', 'Y', 'V']}
        assert _get_motif_profile('N{P}[ST]{P}') == res

    def test__get_motif_profile_empty(self):
        assert _get_motif_profile('') == {}

    