import pytest

from tools.helpers import uni_prot_metadata_parse


class TestTools:

    def test_uni_prot_metadata_parse(self):
        test = '''>sp|P07205|PGK2_HUMAN Phosphoglycerate kinase 2 OS=Homo sapiens OX=9606 GN=PGK2 PE=1 SV=3'''
        assert uni_prot_metadata_parse(test) == {'record_source': 'sp', 'accession': 'P07205', 'entry_name': 'PGK2_HUMAN', 'SV': '3', 'GN': 'PGK2', 'name': 'Phosphoglycerate kinase 2', 'PE': '1', 'OS': 'Homo sapiens', 'OX': '9606'}


    def test_uni_prot_metadata_parse_without_gn(self):
        test = '''>sp|P07205|PGK2_HUMAN Phosphoglycerate kinase 2 OS=Homo sapiens OX=9606 PE=1 SV=3'''
        assert uni_prot_metadata_parse(test) == {'record_source': 'sp', 'accession': 'P07205', 'entry_name': 'PGK2_HUMAN', 'SV': '3', 'GN': None, 'name': 'Phosphoglycerate kinase 2', 'PE': '1', 'OS': 'Homo sapiens', 'OX': '9606'}

    def test_uni_prot_metadata_parse_without_order(self):
        test = '''>sp|P07205|PGK2_HUMAN Phosphoglycerate kinase 2 PE=1 SV=3 OS=Homo sapiens OX=9606'''
        assert uni_prot_metadata_parse(test) == {'record_source': 'sp', 'accession': 'P07205', 'entry_name': 'PGK2_HUMAN', 'SV': '3', 'GN': None, 'name': 'Phosphoglycerate kinase 2', 'PE': '1', 'OS': 'Homo sapiens', 'OX': '9606'}

    def test_uni_prot_metadata_parse_sv_2(self):
        test = '''>sp|P07205|PGK2_HUMAN Phosphoglycerate kinase 2 OS=Homo sapiens OX=9606 GN=PGK2 PE=1 SV=33'''
        assert uni_prot_metadata_parse(test) == {'record_source': 'sp', 'accession': 'P07205', 'entry_name': 'PGK2_HUMAN', 'SV': '33', 'GN': 'PGK2', 'name': 'Phosphoglycerate kinase 2', 'PE': '1', 'OS': 'Homo sapiens', 'OX': '9606'}

    def test_uni_prot_metadata_parse_empty(self):
        with pytest.raises(ValueError):
            uni_prot_metadata_parse('')

    def test_uni_prot_metadata_parse_format(self):
        with pytest.raises(IndexError):
            uni_prot_metadata_parse('ads')