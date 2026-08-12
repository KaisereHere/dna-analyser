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

def convert_to_fasta(fasta_dict):

    '''Converts fasta dictionary to a regular fasta format.

       Args: (dict) fasta_dict: fasta dictionary

       Returns: (str) converted to a fasta format

       Raises: TypeError if fasta_dict has other format than dictionary, ValueError if fasta_dit is empty
    '''
    if not type(fasta_dict) is dict:
        raise TypeError('Fasta has to be a dictionary')
    
    if not fasta_dict:
        raise ValueError("Fasta dictionary can not be empty")

    name = '>' + list(fasta_dict.keys())[0]
    value = list(fasta_dict.values())[0]

    return name + '\n' + str(value)


def uni_prot_metadata_parse(full_name):
    '''Extracts all meta data from the UniProt fasta name

       Args: (str) full_name: fasta name

       Returns: (dict) meta data
    '''

    if not full_name:
        raise ValueError('The name has to be longer than 0')
    
    object_data = {}

    try:
        if full_name.startswith('>'):
            full_name = full_name[1:]
        meta_data = full_name.split()
        source_accession_name = meta_data[0].split('|')
    except IndexError:
        print('The given string does not accord to the UniProt format')
        raise

    object_data['record_source'] = source_accession_name[0]
    object_data['accession'] = source_accession_name[1]
    object_data['entry_name'] = source_accession_name[2]
    object_data['SV'] = None
    object_data['GN'] = None
    object_data['name'] = ''
    object_data['PE'] = None
    object_data["OS"] = None
    object_data["OX"] = None

    rest = ' '.join(meta_data[1:])

    _name = 0
    _values = 1
    _name_end = 2
    mode = _name
    
    key_name = ""
    last_value = ''
    for index, symbol in enumerate(rest):

        if mode == _name:
            object_data['name'] += symbol

            if symbol == "=":     
                mode = _name_end       

        if mode == _values or mode == _name_end: 
            if symbol == "=":
                if key_name:
                    last_value = key_name

                key_name = ''
                key_size = 0
                for inner_symbol in rest[index-1::-1]:

                    if inner_symbol == ' ':
                        break

                    key_size += 1
                    key_name += inner_symbol

                if mode == _name_end:
                    object_data['name'] = object_data['name'][:-key_size-2]
                    mode = _values
                        
                key_name = key_name[::-1]
                object_data[key_name] = ''
                if last_value:
                    object_data[last_value] = object_data[last_value][:-len(key_name)-1]

            else:
                if key_name:
                    object_data[key_name] += symbol

    return object_data