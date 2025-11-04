import functions.ppa_functions as iolo
import pytest
import numpy as np

dict1 = {1: 100, 2: 100, 3: np.nan} 
dict2 = {1: 0, 2: 0, 3: 0, 4: 100}

class TestCombinePdMascotConfs:
    
    def test_combine_conflict_key(self):

        assert iolo.combine_pd_mascot_confs(dict1,dict2)[1] == 100, "The value for key 1 should be 100"

    def test_combine_new_key(self):

        assert iolo.combine_pd_mascot_confs(dict1,dict2)[4] == 100, "The value for the new key should be 100"

    def test_combine_nan(self):

        assert iolo.combine_pd_mascot_confs(dict1,dict2)[3] == 0, "The value for the new key should be 0"

