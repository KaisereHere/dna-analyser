from analyse.sequence import read_fasta, find_components
from tools.helpers import read_file
import graphviz

def main():
    data = read_file('data/rosalind_tree.txt').split('\n')
    
    edges = []
    n_amount = data[0]
    for elem in data[1:]:
        edges.append([int(n) for n in elem.split(' ')])


    ns = [n for n in range(1, int(n_amount)+1)]

    print(len(find_components(ns, edges))-1)

if __name__ == "__main__":
    main()
