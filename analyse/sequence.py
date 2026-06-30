import math


# RNA codon to amino acid table
rna_codon_table = {
    'UUU': 'F','UUC': 'F',
    'UUA': 'L','UUG': 'L','CUU': 'L','CUC': 'L','CUA': 'L','CUG': 'L',
    'AUU': 'I','AUC': 'I','AUA': 'I',
    'AUG': 'M',
    'GUU': 'V','GUC': 'V','GUA': 'V','GUG': 'V',
    'UCU': 'S','UCC': 'S','UCA': 'S','UCG': 'S','AGU': 'S','AGC': 'S',
    'CCU': 'P','CCC': 'P','CCA': 'P','CCG': 'P',
    'ACU': 'T','ACC': 'T','ACA': 'T','ACG': 'T',
    'GCU': 'A','GCC': 'A','GCA': 'A','GCG': 'A',
    'UAU': 'Y','UAC': 'Y',
    'UAA': 'Stop','UAG': 'Stop','UGA': 'Stop',
    'CAU': 'H','CAC': 'H',
    'CAA': 'Q','CAG': 'Q',
    'AAU': 'N','AAC': 'N',
    'AAA': 'K','AAG': 'K',
    'GAU': 'D','GAC': 'D',
    'GAA': 'E','GAG': 'E',
    'UGU': 'C','UGC': 'C',
    'UGG': 'W',
    'CGU': 'R','CGC': 'R','CGA': 'R','CGG': 'R','AGA': 'R','AGG': 'R',
    'GGU': 'G','GGC': 'G','GGA': 'G','GGG': 'G'
}

monoisotopic_mass_table = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333,
    'water': 18.01056
}

def count_nucleotides(strand):

    '''Calculates the amount of the nucleotides in a DNA sequence.

    Args: strand(str): DNA sequence 

    s: Dict with amounts of the nucleotides {A: int, C: int, G: int, T: int} 
    '''
    counted = {'A':0, 'C':0, 'G':0, 'T':0}
    for nucleotide in strand.upper():
        if nucleotide in counted:
            counted[nucleotide] += 1
    return counted

def transcribe(strand):
    '''Simulates transcription of a DNA sequence into a RNA sequence

    Args: strand(str): DNA sequence 

    Returns: (str) RNA sequence 
    '''
    return(strand.upper().replace("T", "U"))

def reversed_complement(strand):
    '''Computes a reversed complementary sequence for a DNA single strand
    
    Args: strand(str): DNA sequence 

    Returns: (str) DNA sequence 
     
    '''
    return strand[::-1].upper().replace("A", "t").replace("T", "a").replace("G","c").replace("C", "g").upper() 

def read_fasta(fasta):
    '''Converts the given fasta content into a python dict 

    Args: fasta: fasta content >name
                               GGACT

    Returns: dict {'name': 'GGACT'} with the named sequences imported from a fasta content  
    '''
    strands_dataset = {}

    current_label = ''

    for item in fasta.splitlines():
        if item:
            if item.startswith('>'):
                current_label = item[1::]
                strands_dataset[current_label] = ''
        
            else: strands_dataset[current_label] += item
    
    return strands_dataset

def calculate_gc_content(strand):

    '''Calculates the percentage of G and C nucleotides in a DNA sequence

        Args: strand(str): DNA sequence 

        Returns: (float) Percentage 
    '''
    all_nucleotides = count_nucleotides(strand)

    total = (
        all_nucleotides['A'] +
        all_nucleotides['C'] +
        all_nucleotides['G'] +
        all_nucleotides['T']
    )

    return (all_nucleotides['G'] + all_nucleotides['C']) / total * 100 if total else 0

def calculate_hamming_distance(strand1, strand2):
    '''Calculates hamming distance between two DNA sequences 

        Args: strand1(str): DNA sequence , strand2(str): another DNA sequence 

        Returns: (int) distance 

        Raise: ValueError if sequences have not the same length 
    '''
    distance = 0

    if len(strand1) != len(strand2):
      raise ValueError(f"Sequences must be equal length: {len(strand1)} vs {len(strand2)}")
    
    for nucleotide1, nucleotide2 in zip(strand1, strand2):
        if nucleotide1 != nucleotide2:
            distance += 1

    return distance

def find_motif(strand, motif):
    '''Finds all places where the given motif occurs

        Args: strand(str): DNA sequence , motif(str): another DNA sequence 

        Returns: list with indices (starts with 1 - string [1] = s)  
    '''
    motifs = []
    if motif != '' and strand != '':
        for nucleotide_number in range(len(strand)):
            if strand[nucleotide_number:len(motif) + nucleotide_number] == motif:
                motifs.append(nucleotide_number+1)
    return motifs

def translate(strand):
    '''Simulates translation of a RNA sequence into an amino acid sequence (protein)

        Args: strand(str): RNA sequence 

        Returns: (str) amino acid sequence (protein) 

        Raise: ValueError if the unknown codon is found 
    '''
    protein = []
    strand = strand.strip().upper()

    for codon_number in range(0, len(strand), 3):
        codon = strand[codon_number:codon_number+3]

        if codon not in rna_codon_table or len(codon) < 3:
            raise ValueError('The given codon does not exist')

        if rna_codon_table[codon] == 'Stop':
            break
        
        protein.append(rna_codon_table[codon])

    return ''.join(protein)

    
def gc_sliding_window(sequence, window_size):

    '''Examines the local GC contents on the subsequence of the given DNA or RNA sequence 
    
        Args: sequence(str): DNA or RNA sequence, window_size(int): the size of the subsequences 

        Returns: list with GC percentages of GC content

        Raise: IndexError if the window size is not appropriate 
    '''
    gc_list = []

    if window_size > len(sequence) or window_size <= 0:
        raise IndexError('The window size is larger than a length of the sequence or equal 0')
    
    for nucleotide_number in range(0, len(sequence)-window_size+1):
        gc_list.append(calculate_gc_content(sequence[nucleotide_number:nucleotide_number+window_size]))

    return gc_list

def dominant_probability(dominant_homozygous, heterozygous, recessive_homozygous):

    '''Calculates the probability that an offspring exhibits the dominant phenotype

        Args: dominant_homozygous(int):  amount of dominant homozygous 
              heterozygous(int):         amount of heterozygous 
              recessive homozygous(int): amount of recessive homozygous 

        Returns: (float) probability

        Raise: ValueError if one of the given arguments is negative
               TypeError  if one of the given arguments is not an integer
    '''
    if dominant_homozygous < 0 or heterozygous < 0 or recessive_homozygous < 0:
        raise ValueError("A size of a population can not be negative")
    
    if type(dominant_homozygous) is not int or type(heterozygous) is not int or type(recessive_homozygous) is not int:
        raise TypeError("The amount of organisms supossed to be an integer")
    
    sum_organisms = dominant_homozygous + heterozygous + recessive_homozygous

    if sum_organisms < 2:
        return 0

    two_hetero = heterozygous/sum_organisms * (heterozygous-1)/(sum_organisms-1) * 1/4
    two_recessive = recessive_homozygous/sum_organisms * (recessive_homozygous-1)/(sum_organisms-1)
    heterozygous_recessive = heterozygous/sum_organisms * recessive_homozygous/(sum_organisms-1) * 2/4
    recessive_hetero = recessive_homozygous/sum_organisms * heterozygous/(sum_organisms-1) * 2/4
    
    return 1 - two_hetero - two_recessive - heterozygous_recessive - recessive_hetero
    
    
def find_monoisotopic_mass(sequence):
    """Calculates monoisotopic mass of the given protein 

        Args: sequence(str): protein 

        Returns: (float) mass 

        Raise: ValueError if an unknown amino acid found
    """
    mass = 0

    for amino_acid in sequence:

        if amino_acid not in monoisotopic_mass_table:
            raise(ValueError(f'"{amino_acid}" - such amino acid name is not assigned'))
        
        mass += monoisotopic_mass_table[amino_acid]

    return mass


def strand_profile(strands):

    """Builds profile of DNA strand based on comparing its differently mutated versions. Counts the occurrence of each nucleotide
        on the corresponding location in the given sequences.

        Args: strands(custom fasta dict): different versions of DNA 

        Returns: (custom dict) profile 

        Raise: ValueError if list is empty, if strands have not the same length or if a nucleotide is found which does not exist 
    """

    if not strands:
        raise ValueError("The given list can not be empty")

    strand_size = len(next(iter(strands.values())))
    
    profile = {'A': [0] * strand_size,'G': [0] * strand_size, 'C': [0] * strand_size,'T': [0] * strand_size}

    for strand_name, sequence in strands.items():
        
        if len(sequence) != strand_size:
            raise ValueError(f"The length of the {strand_name} strand is not equal to other")
        
        for index, nucleotide in enumerate(sequence):

            if nucleotide not in profile:
                raise ValueError(f'The nucleotide "{nucleotide}" does not exist')
                       
            profile[nucleotide][index] += 1

    return profile


def consensus_strand(profile):
    """Finds out the most probable DNA sequence (consenus) by analysing it's profile

        Args: profile(custom dict): profile of a DNA strand 

        Returns: (str)consensus 
    """
    strand_length = len(profile['A']) 

    if strand_length == 0:
        return ''
    
    strand = ['a'] * strand_length

    for amount_index in range(strand_length):
        strand[amount_index] = max(profile, key=lambda nuc: profile[nuc][amount_index])

    return ''.join(strand)
 
def open_read_frame(original_strand):
    
    '''Finds all possible proteins given on both strands DNA sequence

        Args: original_strand(str): DNA strand 

        Returns: (set) all possible amino acid sequences 
    '''
    proteins = []
    complement_strand = transcribe(reversed_complement(original_strand))
    original_strand = transcribe(original_strand)

    for strand in (original_strand, complement_strand):
        for shift in range(3):

            for codon_number in range(shift, len(strand)-(len(strand)-shift)%3, 3): # disignate the reading after shifting

                triplet = strand[codon_number:codon_number+3]
                if triplet == 'AUG':
    
                    for translation_codon_number in range(codon_number,  len(strand)-(len(strand)-shift)%3, 3): # begin of the protein  

                        triplet_translation = strand[translation_codon_number:translation_codon_number+3] # protein triplet 
                        if rna_codon_table[triplet_translation] == 'Stop':
                            proteins.append(translate(strand[codon_number:translation_codon_number+3])) # complete protein
                            break
            

    return set(proteins)


def rna_splicing(rna_strand, introns):

    '''Removes introns from RNA strand

        Args: rna_strand(str): RNA sequence , introns(list with str): list with introns 

        Returns: (str) spliced RNA sequence 

        Raise: ValueError if intron is longer than a strand
    '''
    for intron in introns:
        if len(intron) > len(rna_strand):
            raise ValueError('Intron can not be larger than a strand')
        rna_strand = rna_strand.replace(intron, '')

    return rna_strand

def dominant_offspring(AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa):
    """Calculates an expected amount of dominant offsprings based on the given pairs

        Args: int: AA_AA: amount of dominant homozygous pairs, AA_Aa: dominant homozygous and heterozygous, Aa_Aa: heterozygous
                   Aa_aa: heterozygous and recessive homozygous, aa_aa: recessive homozygous

        Returns: (int) Amount of offsprings 
    """
    AA_AA_assumption = 2
    AA_Aa_assumption = 2
    AA_aa_assumption = 2
    Aa_Aa_assumption = 6/4
    Aa_aa_assumption = 1
    aa_aa_assumption = 0

    return (AA_AA_assumption*AA_AA + AA_Aa_assumption*AA_Aa + AA_aa_assumption*AA_aa +
            Aa_Aa_assumption*Aa_Aa + Aa_aa_assumption*Aa_aa + aa_aa_assumption*aa_aa)

def is_substring(sequence_dict, substring):
    """Checks whether the every DNA sequence in a dict consists the given substring

        Args: sequence_dict(custom fasta dict): DNA sequences, substring(str): substring

        Returns: (bool) True or False 
    """
    for sequence in sequence_dict.values():
        if substring not in sequence:
            return False

    return True    

def find_shared_motif(sequence_dict): 
    '''Finds a motif which occurs in every DNA sequence

        Args: sequence_dict(custom fasta dict): DNA sequences 

        Returns: (str) motif 

        Raise: ValueError if the given dictionary is empty
    '''
    if len(sequence_dict) == 0:
        raise ValueError("The dictionary can not be empty")
    
    shortest_item = min(sequence_dict.values(), key=lambda x:len(x))

    for shift in range(len(shortest_item), -1, -1):
        for index in range(len(shortest_item)-shift):

            substring = shortest_item[index:index + shift + 1]

            if is_substring(sequence_dict, substring):
                return substring
            
    return ''
            
def independent_alleles(generation, count_heterozygous):
    '''Calculates the birth probability of a dominantly phenotypical expressed offspring in the n-generation in a k-size population
        (Mendelian segregation law) 

        Args: generation(int): number of generetion, count_heterozygous(int): size of a heterozygous popuation 

        Returns: (float) probability 
    '''
    tries = pow(2, generation)
    probability = 0

    for count in range(count_heterozygous, tries+1):
        probability += (math.factorial(tries) / (math.factorial(count) * math.factorial(tries-count)))  * ((pow(1/4, count) * pow(3/4, tries-count)))

    return probability

def spliced_motif(sequence, subsequence):

    '''Checks whether the given sequence consists the subsequence after deletion of introns.

        Args: sequence(str): DNA sequence, subsequence(str):another DNA sequence  

        Returns: (list) of occurrence indices 

        Raise: ValueError if the given subsequence is empty
    '''
    if len(subsequence) == 0:
        raise ValueError("The subsequence can not be empty")

    indices = []

    subsequence_index = 0

    for index in range(len(sequence)):
        
        
        if sequence[index] == subsequence[subsequence_index]:
            indices.append(index)
            subsequence_index += 1

            if len(subsequence) == subsequence_index:
                return indices

    return []
    

def shared_spliced_motif(sequence1, sequence2):

    """Uses LCS to find the longest subsequence in the 2 DNA sequences 
        Complexity: O(n×m)

        Args: sequence1(str): DNA sequence, sequence2(str):another DNA sequence 

        Returns: (str) subsequence 

    """
    longest_subsequences = [[0 for _ in range(len(sequence2) + 1)] for _ in range(len(sequence1) + 1)]
    res = ''
    for index1 in range(len(sequence1)-1, -1, -1):
        for index2 in range(len(sequence2)-1, -1, -1):

            if sequence1[index1] == sequence2[index2]:
                longest_subsequences[index1][index2] = longest_subsequences[index1+1][index2+1] + 1
            else:
                longest_subsequences[index1][index2] = max(longest_subsequences[index1+1][index2], longest_subsequences[index1][index2+1])

    i, j = 0, 0

    while i < len(sequence1) and j < len(sequence2):
        if sequence1[i] == sequence2[j]:
            res += sequence1[i]
            i += 1
            j += 1

        else:
            if longest_subsequences[i+1][j] > longest_subsequences[i][j+1]:
                i += 1
            else:
                j += 1


    return res
    
def edit_distance(sequence1, sequence2):

    """Uses Levenshtein's algorithm to find the minimum number of single-character edits required to transform one string into another.
        Complexity: O(n×m)

        Args: sequence1(str): DNA sequence, sequence2(str):another DNA sequence 

        Returns: (int) the minimum number 
    """
    shortest_path = [[0 for _ in range(len(sequence2) + 1)] for _ in range(len(sequence1) + 1)]

    length1 = len(sequence1)
    length2 = len(sequence2)

    for index1 in range(length1):
        shortest_path[index1][length2] = length1 - index1
    
    for index2 in range(length2):
        shortest_path[length1][index2] = length2 - index2

    for index1 in range(length1-1, -1, -1):
        for index2 in range(length2-1, -1, -1):

            if sequence1[index1] == sequence2[index2]:
                shortest_path[index1][index2] = shortest_path[index1+1][index2+1]

            else:
                shortest_path[index1][index2] = min(shortest_path[index1+1][index2], shortest_path[index1][index2+1], shortest_path[index1+1][index2+1]) + 1 
        
    return shortest_path[0][0] 