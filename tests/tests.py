"""

Bioinformatics BIOL68400 Programming Assignment: Tests for LRG (XML) file parser
Authors: Callum Rakhit and Graeme Smith
Created: 21st November 2018
Description: Contains tests for lrg_parser and associate scripts
Usage: See README for detailed documentation.

"""

import unittest2 as unittest

from LRG_parser import get_summary_list, get_lrg_id

class lrgParserTest(unittest.TestCase):
    """Tests for LRG_parser.py."""

    def test_get_lrg_id(self):
        """Is LRG ID returned correctly"""
        self.assertEquals(get_lrg_id("COL1A1", "hgnc"),"LRG_1",
                          msg="ERROR: ID returned by get_lrg_id() does not match expected output")
        self.assertEquals(get_lrg_id("COL1A2", "hgnc"),"LRG_2",
                          msg="ERROR: ID returned by get_lrg_id() does not match expected output")
        self.assertEquals(get_lrg_id("NG_007400.1", "xref"),"LRG_1",
                          msg="ERROR: ID returned by get_lrg_id() does not match expected output")
        self.assertEquals(get_lrg_id("NG_007405.1", "xref"),"LRG_2",
                          msg="ERROR: ID returned by get_lrg_id() does not match expected output")

# TODO write relevant tests

if __name__ == '__main__':
    unittest.main()

#
# def test_get_summary(tree):
#     """
#     Checks the output of get_summary() against expected results.
#     This test assumes you picked LRG_1 so only works for one use case.
#     """
#     expected_results = [('Gene', 'COL1A1'), ('LRG version', 'LRG_1'),
#                         ('HGNC ID', '2197'), ('Sequence source', 'NG_007400.1')]
#     assert get_summary_list(tree) == expected_results, 'Summary results not as expected!'
#     print("get_summary(): passed")

