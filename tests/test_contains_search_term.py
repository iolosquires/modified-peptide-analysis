import functions.ppa_functions as iolo
import pytest
acc_list = ['1123','456','789']


@pytest.mark.parametrize("pass_search", ['1123','456'])

def test_in_list(pass_search):
    assert iolo.contains_search_term(acc_list,pass_search) == True, "Expected to be found in the list"

@pytest.mark.parametrize("fail_search", ['123','56','78','112'])

def test_not_in_list(fail_search):
    assert iolo.contains_search_term(acc_list,fail_search) == False, "Expected not to be found in the list"

