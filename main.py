from analyse.sequence import read_fasta, find_motif_regex, find_motif, _get_motif_profile
from tools.helpers import read_file, get_protein_uniprot, parse_template_salt
import time

def main():
    res = _get_motif_profile('N{P}[ST]{P}')
    print(res)


if __name__ ==  '__main__':
    main()