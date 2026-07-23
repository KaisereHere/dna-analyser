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

def _uni_prot_metadata_parse(full_name):
    '''Extracts all meta data from the UniProt fasta name

       Args: (str) full_name: fasta name

       Returns: (dict) meta data
    '''
    object_data = {}

    meta_data = full_name.split()
    source_accession_name = meta_data[0].split('|')

    object_data['record_source'] = source_accession_name[0]
    object_data['accession'] = source_accession_name[1]
    object_data['entry_name'] = source_accession_name[2]

    rest = ' '.join(meta_data[1:])

    mode = -1
    object_data['ox'] = ''
    object_data['gn'] = ''
    object_data['name'] = ''
    object_data['os'] = ''

    for index, symbol in enumerate(rest):

        if mode == -1:
            object_data['name'] += symbol


        if mode == 1:
            
            if symbol == " ":
                mode = 0
                continue

            object_data['ox'] += symbol

        if mode == 2:
            if symbol == " ":
                mode = 0
                continue

            object_data['gn'] += symbol

        if mode == 3:
            object_data['os'] += symbol
                
        if symbol == '=':
            service = rest[index-2:index]

            if service == 'OS':
                mode = 3
                object_data['name'] = object_data['name'][:index-3]

            if service == 'OX':
                mode = 1
                object_data['os'] = object_data['os'][:-4]

            if service == 'GN':
                mode = 2

            if service == "PE":
                object_data['pv'] = rest[index+1] 

            if service == "SV":
                object_data['sv'] = rest[index+1] 

    return object_data