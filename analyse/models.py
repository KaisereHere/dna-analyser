from analyse.sequence import read_fasta, translate, reversed_complement as rc 
from analyse.sequence import calculate_gc_content, count_nucleotides, transcribe
from analyse.sequence import find_motif_regex, _get_motif_profile, open_read_frame
from analyse.sequence import find_monoisotopic_mass, rna_splicing

from tools.helpers import _uni_prot_metadata_parse
from dataclasses import dataclass

@dataclass
class Sequence:
    sequence: str

@dataclass
class FastaSequence:
    fasta: str

    def __post_init__(self):
        self.fasta_dict = read_fasta(self.fasta)

    @property
    def sequence(self):
        return list(self.fasta_dict.values())[0]

    @property
    def full_name(self):
        return list(self.fasta_dict.keys())[0]


@dataclass
class NucleicAcidSequence(Sequence):
    @property
    def nucleotide_count(self):
        return count_nucleotides(self.sequence)

@dataclass
class RNASequence(NucleicAcidSequence):

    @property
    def translation(self):
        return ProteinSequence(translate(self.sequence))

    def splice(self):
        return RNASequence(rna_splicing(self.sequence))


@dataclass 
class DNASequence(NucleicAcidSequence):

    @property
    def transcription(self):
        return RNASequence(transcribe(self.sequence))
    
    @property
    def translation(self):
        return self.transcription.translation

    @property
    def reverse_complement(self):
        return DNASequence(rc(self.sequence))

    @property
    def gc_content(self):
        return calculate_gc_content(self.sequence)

    def find_motif(self, motif):
        return self.translation.find_motif(motif)

    @property
    def all_proteins(self):
        return [ProteinSequence(sequence) for sequence in open_read_frame(self.sequence)]


@dataclass
class ProteinSequence(Sequence):

    @property
    def monoisotopic_mass(self):
        return find_monoisotopic_mass(self.sequence)

    def find_motif(self, motif):
        profile = _get_motif_profile(motif)
        return find_motif_regex(profile, self.sequence)

@dataclass
class UniProtRecord:

    fasta: FastaSequence

    
    def __post_init__(self):
        object_data = _uni_prot_metadata_parse(self.fasta.full_name)

        self.sequence = ProteinSequence(self.fasta.sequence)

        self.accession = object_data['accession']
        self.entry_name = object_data['entry_name']
        self.taxonomy_id = object_data['ox']
        self.gene_name = object_data['gn']
        self.organism = object_data['os']
        self.sequence_version = object_data['sv']
        self.protein_existance = object_data['pv']
        self.protein_name = object_data['name']
