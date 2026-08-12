from analyse.sequence import read_fasta, translate, reversed_complement as rc 
from analyse.sequence import calculate_gc_content, count_nucleotides, transcribe
from analyse.sequence import find_motif_regex, prosite_parse, open_read_frame
from analyse.sequence import find_monoisotopic_mass, rna_splicing, edit_distance
from analyse.sequence import strand_profile, consensus_strand, calculate_hamming_distance
from analyse.sequence import align as _align

from tools.helpers import uni_prot_metadata_parse, convert_to_fasta
from dataclasses import dataclass

@dataclass
class Sequence:
    sequence: str

    def align(self, seq2:Sequence|str):
        
        if type(seq2) == str:
            seq2 = Sequence(seq2)

        alignment1, alignment2, _ = _align(self.sequence, seq2.sequence)
        return (type(self)(alignment1), type(self)(alignment2))

    def edit_distance(self, seq2:Sequence|str):
        
        if type(seq2) == str:
            seq2 = Sequence(seq2)

        return edit_distance(self.sequence, seq2.sequence)

    def hamming_distance(self, seq2:Sequence|str):
        if type(seq2) == str:
            seq2 = Sequence(seq2)

        return calculate_hamming_distance(self.sequence, seq2.sequence)
    
    @classmethod
    def consensus(cls, sequences):

        if type(sequences) == list:
            sequences = dict(enumerate(sequences))

        if isinstance(list(sequences.values())[0], Sequence):
            sequences = {name: value.sequence for name, value in sequences.items()}


        consensus_profile = strand_profile(sequences)
        consensus_seq = consensus_strand(consensus_profile)

        return cls(consensus_seq)

@dataclass
class FastaSequence:
    fasta: str | dict

    @classmethod
    def array(cls, fasta:str | dict):

        array = []
        if isinstance(fasta, str):
            fasta = read_fasta(fasta)

        for name, sequence in fasta.items():
            array.append(cls({name:sequence}))

        return array
    
    def __post_init__(self):
        arg_type = type(self.fasta)

        if arg_type is dict:
            self.fasta = convert_to_fasta(self.fasta)

        elif arg_type is FastaSequence:
            self.fasta = self.fasta.fasta

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

    def splice(self, introns):
        return RNASequence(rna_splicing(self.sequence, introns))

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

    def find_protein_motif(self, motif):
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
        profiles = prosite_parse(motif)
        matches = []
        for profile in profiles:
            matches.extend(find_motif_regex(profile, self.sequence))
        return matches

    
class UniProtRecord:

    def __init__(self, fasta=str|dict|FastaSequence):
        self.fasta = FastaSequence(fasta)
        object_data = uni_prot_metadata_parse(self.fasta.full_name)
        self.sequence = ProteinSequence(self.fasta.sequence)

        self.accession = object_data['accession']
        self.entry_name = object_data['entry_name']
        self.taxonomy_id = object_data['OX']
        self.gene_name = object_data['GN']
        self.organism = object_data['OS']
        self.sequence_version = object_data['SV']
        self.protein_existence = object_data['PE']
        self.protein_name = object_data['name']


