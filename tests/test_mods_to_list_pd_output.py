
import functions.ppa_functions as iolo
pd_string = 'S20(Phospho): 99.97'

class TestClass:
    def test_one(self):
        assert iolo.mods_to_list_pd_output(pd_string) == [['99.97', 'S20', 'Phospho)']]



