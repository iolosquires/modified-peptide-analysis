import functions.ppa_functions as iolo
import pytest

@pytest.mark.parametrize("single", ['PEPTIDE','PEPYIDE','PEPSIDE'])

def test_single(single):
    assert iolo.find_phospho_single_acceptor_site(single) == True

@pytest.mark.parametrize("not_single", ['PEPTSIDE','PEYYIDE','PEPSTIDE'])

def test_not_single(not_single):
    assert iolo.find_phospho_single_acceptor_site(not_single) == False


    