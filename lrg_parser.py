"""
Authors:           Callum Rakhit and Graeme Smith
Created:           November 2018
Description:       Parses LRG XML file and outputs BED file
Usage:             See README for detailed documentation.
"""

# Python 3.5
import xml.etree.ElementTree as ET
#import pandas as pd
import warnings
import os
import glob
import sys

# Download relevant LRG file

# Hard coded test file:
lrg_file = "/home/graeme/Desktop/UoM_LRG_Parser/LRG_public_xml_files/LRG_1.xml"

# Import data from LRG file

tree = ET.parse(lrg_file)

root = tree.getroot()

for child in root:
    print(child.tag, child.attrib)

print([elem.tag for elem in root.iter()])

# Write imported data to BED file

def lrg_2_bed():
    """Converts """