import requests

amino_acids_one_letter = [
    'A', 'R', 'N', 'D', 'C', 
    'Q', 'E', 'G', 'H', 'I', 
    'L', 'K', 'M', 'F', 'P', 
    'S', 'T', 'W', 'Y', 'V'
]

def read_file(filename):
    with open(filename) as f:
        return f.read()
    
def get_protein_uniprot(id):
    '''Provides the taken from uniprot.org protein sequence in fasta format.
    
       Args: id(str): protein id for the uniprot query

       Returns: (str)raw fasta text
    '''
    id = id.split('_')[0]
    res = requests.get(f'https://rest.uniprot.org/uniprotkb/{id}.fasta')
    return res.text

def parse_template_salt(template):
    _normal = 1
    _or = 2
    _except = 3

    _modus = _normal
    
    profile = {}

    counter = 0
    excepts = []

    for symbol in template:
        if symbol == '{':
            _modus = _except
            continue
            
        if symbol == '}':
            _modus = _normal
            for amino_acid in amino_acids_one_letter:
                if amino_acid not in excepts:
                    profile[counter].append(amino_acid)

            excepts = []
            counter += 1
            continue

        if symbol == '[':
            _modus = _or
            continue

        if symbol == ']':
            _modus = _normal
            counter += 1
            continue

        if _modus == _normal:
            profile[counter] = [symbol]
            counter += 1
            continue

        if _modus == _or:
            if counter not in profile:
                profile[counter] = []
            
            profile[counter].append(symbol)
            continue

        if _modus == _except:
            if counter not in profile:
                profile[counter] = []

            excepts.append(symbol)

    variants = ['']
    new = []

    for position, lst in profile.items():
        for variant in variants:
            for symbol in lst:
                new.append(variant+symbol)
        variants = new
        new = []
            
    return variants