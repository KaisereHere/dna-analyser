
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