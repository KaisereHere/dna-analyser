from analyse.models import UniProtSequence
from tools.helpers import get_protein_uniprot

def main():
    data = get_protein_uniprot('P07204')
    test = UniProtSequence(data)
    print(test)

if __name__ == "__main__":
    main()