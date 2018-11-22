"""

Bioinformatics BIOL68400 Programming Assignment: Tests for LRG (XML) file parser
Authors: Callum Rakhit and Graeme Smith
Created: 21st November 2018
Description: Contains tests for lrg_parser and associate scripts
Usage: See README for detailed documentation.

"""

# import unittest2 as unittest
# from helpers import *
from UoM_LRG_Parser.LRG_parser import get_summary_list

# TODO write relevant tests

if __name__ == '__main__':
    unittest.main()


def test_get_summary(tree):
    """
    Checks the output of get_summary() against expected results.
    This test assumes you picked LRG_1 so only works for one use case.
    """
    expected_results = [('Gene', 'COL1A1'), ('LRG version', 'LRG_1'),
                        ('HGNC ID', '2197'), ('Sequence source', 'NG_007400.1')]
    assert get_summary_list(tree) == expected_results, 'Summary results not as expected!'
    print("get_summary(): passed")
