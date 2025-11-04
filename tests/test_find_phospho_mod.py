import functions.ppa_functions as iolo

mod_dict_list = [{'location': 3, 'name': 'Phospho'},{'location': 4, 'name': 'Phospho'},{'location': 4, 'name': 'Oxidation'}]

class TestClass:
    def test_one(self):
        assert iolo.find_phospho_mod(mod_dict_list) == True

    