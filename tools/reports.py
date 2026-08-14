from analyse.models import UniProtRecord
from tools.helpers import get_protein_uniprot

def protein_report(protein):

    protein = UniProtRecord(get_protein_uniprot(protein))
    print("=====================")
    print("Protein summary")
    print("=====================")

    print('Uniprot metadata:\n')
    print('Accession:', protein.accession)
    print('Entry name:', protein.entry_name)
    print("Gene name:", protein.gene_name)
    print("Protein name:", protein.protein_name)
    print("Organism:", protein.organism)
    print("Sequence version:", protein.sequence_version)
    print("Protein existence:", protein.protein_existence)

    print("=====================")
    print("General information:\n")    
    print("Protein sequence:", protein.sequence.sequence[:10] + "... total length:", len(protein.sequence.sequence))
    print("Monoisotopic mass:", protein.sequence.monoisotopic_mass)
    print("=====================")
    print('Motifs:\n')
    glycosilation_sites = protein.sequence.find_motif(r'N-{P}-[ST]-{P}')
    print("Potential N-glycosilation sites:", glycosilation_sites, "overall amount: ", len(glycosilation_sites))

    phosphorylation_sites = protein.sequence.find_motif('[ST]-x(2)-[DE]')
    print("Potential N-phosphorylation sites:", phosphorylation_sites, "overall amount: ", len(phosphorylation_sites))
    print("=====================\n")


def compare_2_proteins_report(protein1, protein2):
    print("=====================")
    print("Protein comparation")
    print("=====================")
    protein1 = UniProtRecord(get_protein_uniprot(protein1))
    protein2 = UniProtRecord(get_protein_uniprot(protein2))

    res = protein1.sequence.align(protein2.sequence)

    print(f'{protein1.protein_name}\t {res[0].sequence}\n{protein2.protein_name} {res[1].sequence}')
    edit_distance = protein1.sequence.edit_distance(protein2.sequence)
    print(f'Edit distance:', edit_distance, "amino acids")

    identity = len(res[0].sequence) - res[0].hamming_distance(res[1])  

    print("Identity:", identity, "amino acids,", str(identity/len(res[0].sequence)*100), "%")
    print("=====================")
