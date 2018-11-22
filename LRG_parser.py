"""
Authors:           Callum Rakhit and Graeme Smith
Created:           November 2018
Description:       Parses LRG XML files and outputs a BED file
Usage:             See README for detailed documentation.
"""

# Python 3.5
import argparse
import glob
import os
import pandas as pd
import requests
import sys
import warnings
import xml.etree.ElementTree as ET


""" 
Import Arguments from command line
"""
parser = argparse.ArgumentParser(
    description='Downloads and parses Locus Reference Genomic (LRG) files and produces a BED file')
parser.add_argument('-hgnc', '-h',
                    nargs='+',
                    help='Import LRG files for conversion into a BED file as per provided HGNC IDs',
                    required=False)
parser.add_argument('-xref', '-x',
                    nargs='+',
                    help='Import LRG files for conversion into a BED file as per provided external references',
                    required=False)
parser.add_argument('-output_dir', '-o',
                    type=str,
                    help='File path to output BED file. Defaults to current working directory',
                    required=False)
parser.add_argument('-bed_file', '-b',
                    type=str,
                    help='Specify the name of the BED file which the script will output. '
                         'If not specified will automatically append .bed file suffix',
                     required=False)
args = parser.parse_args()


def get_valid_lrg_id_list():
    """Uses the EMBL-EBI EB-eye RESTful service to retrieve list of valid files"""
    pass

def fetch_lrg_file(ref, type):
    pass


def get_lrg_id(ref, type):
    """Uses the EMBL-EBI RESTful API to translate hgnc symbols or external refs into LRG file name"""
    if type == "hgnc":
        url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=name:' + ref
    elif type == "xref":
        url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=' + ref
    response = requests.get(url, allow_redirects=True)
    root = ET.fromstring(response.content)
    # Parse the returned xml file for the LRG file name
    for entry in root.iter('entry'):
        lrg_id = entry.attrib["id"]
    return lrg_id


def lrg2bed():
    """Parses LRG format files and converts to BED file format"""
    pass




def import_lrg_files():
    """Imports multiple requested LRG files"""
    pass


def check_lrg_id(lrg_id):
    """Checks that an lrg_id is valid"""
    pass



def write_bed_file():
    """Writes a BED file"""
    pass





