
import functions.ppa_functions as iolo

class TestClass:
    def test_one(self):
        assert iolo.get_phospho_position('030') == [0]

    def test_two(self):
        assert iolo.get_phospho_position('0340') == [0,1]

    def test_three(self):
        assert iolo.get_phospho_position('040') == [0]

