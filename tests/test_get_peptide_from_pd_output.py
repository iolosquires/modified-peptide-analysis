import functions.ppa_functions as iolo
import pytest

@pytest.mark.parametrize("single", ['ASDsadfs.PEPTIDE.asdsadsad'])

def test_single(single):
    assert iolo.get_peptide_from_pd_output(single) == 'PEPTIDE'

