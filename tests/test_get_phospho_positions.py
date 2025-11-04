import functions.ppa_functions as iolo

mod_dict_list = [{'location': 3, 'name': 'Phospho'},{'location': 4, 'name': 'Phospho'},{'location': 4, 'name': 'Oxidation'}]

class TestClass:
    def test_one(self):
        assert iolo.get_phospho_positions(mod_dict_list) == '3p4p', "Expected 3p4p"

    