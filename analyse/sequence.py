from re import finditer

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
    counted = {'A':0, 'C':0, 'G':0, 'T':0}
    for nucleotide in strand.upper():
        if nucleotide in counted:
            counted[nucleotide] += 1
    return counted

def transcribe(strand):
    return(strand.upper().replace("T", "U"))

def reversed_complement(strand):
    return strand[::-1].upper().replace("A", "t").replace("T", "a").replace("G","c").replace("C", "g").upper() 

def read_fasta(fasta):

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

    all_nucleotides = count_nucleotides(strand)

    total = (
        all_nucleotides['A'] +
        all_nucleotides['C'] +
        all_nucleotides['G'] +
        all_nucleotides['T']
    )

    return (all_nucleotides['G'] + all_nucleotides['C']) / total * 100 if total else 0

def calculate_hamming_distance(strand1, strand2):

    distance = 0

    if len(strand1) != len(strand2):
      raise ValueError(f"Sequences must be equal length: {len(strand1)} vs {len(strand2)}")
    
    for nucleotide1, nucleotide2 in zip(strand1, strand2):
        if nucleotide1 != nucleotide2:
            distance += 1

    return distance

def find_motif(strand, motif):
    motifs = []
    if motif != '' and strand != '':
        for nucleotide_number in range(len(strand)):
            if strand[nucleotide_number:len(motif) + nucleotide_number] == motif:
                motifs.append(nucleotide_number+1)
    return motifs

def translate(strand):
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

    gc_list = []

    if window_size > len(sequence) or window_size <= 0:
        raise IndexError('The window size is larger than a length of the sequence or equal 0')
    
    for nucleotide_number in range(0, len(sequence)-window_size+1):
        gc_list.append(calculate_gc_content(sequence[nucleotide_number:nucleotide_number+window_size]))

    return gc_list

def dominant_probability(dominant_homozygous, heterozygous, recessive_homozygous):

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

    mass = 0

    for amino_acid in sequence:

        if amino_acid not in monoisotopic_mass_table:
            raise(ValueError(f'"{amino_acid}" - such protein name is not assigned'))
        
        mass += monoisotopic_mass_table[amino_acid]

    return mass


def strand_profile(strands):

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
    strand_length = len(profile['A']) 

    if strand_length == 0:
        return ''
    
    strand = ['a'] * strand_length

    for amount_index in range(strand_length):
        strand[amount_index] = max(profile, key=lambda nuc: profile[nuc][amount_index])

    return ''.join(strand)
 
def open_read_frame(original_strand):
    
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