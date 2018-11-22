"""

Bioinformatics BIOL68400 Programming Assignment: Parser for LRG (XML) files
Authors: Callum Rakhit and Graeme Smith
Created: 21st November 2018
Description: Parses LRG XML files and outputs a BED file
Usage: See README and documentation.docx for detailed documentation.

"""

# Python 3.5/3.6.6

from urllib.request import urlopen  # Python 3 specific
import argparse
# import glob
import os
# import pandas as pd
import requests
import sys
# import warnings
import xml.etree.ElementTree as ET


""" 
Import Arguments from command line
"""

parser = argparse.ArgumentParser(
    description='Downloads and parses Locus Reference Genomic (LRG) files and produces a BED file')
file_location = parser.add_mutually_exclusive_group(required=True)
file_location.add_argument('-l', '--local', type=int, action='store',
                    help='Takes the LRG ID and parses a copy of the lrg'
                         'file from the local directory. Assumes file is using'
                         'the same naming convention as the LRG website.'
                         'i.e. LRG_{user input}.xml')
file_location.add_argument('-w', '--web', type=int, action='store',
                               help='Takes the LRG ID and parses a copy of the LRG file '
                               'from the LRG website')
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


def get_args():
    """
    Function uses argparse to setup command line arguments and catch errors
    """
    parser = argparse.ArgumentParser()
    file_location = parser.add_mutually_exclusive_group(required=True)
    file_location.add_argument('-l', '--local', type=int, action='store',
                               help='Takes the LRG ID and parses a copy of the lrg'
                                    'file from the local directory. Assumes file is using'
                                    'the same naming convention as the LRG website.'
                                    'i.e. LRG_{user input}.xml')
    file_location.add_argument('-w', '--web', type=int, action='store',
                               help='Takes the LRG ID and parses a copy of the LRG file '
                               'from the LRG website')
    # parser.add_argument('-f', '--file', action='store_true',
    #                     help="Optional: writes the output to a file instead of to the console.")
    return parser.parse_args()


def get_lrg_file(sys_args):
    """
    This pulls the LRG/XML file either locally or from the LRG FTP site and
    parses using ElementTree. sys_args.local means local file and sys_args.web
    means retrieve file from the LRG website.
    """
    if sys_args.local:  # If the user specified a local file
        try:
            lrg_xml = open('LRG_' + str(sys_args.local) + '.xml', 'r')  # Use the local file
            print('LRG_' + str(sys_args.local) + '.xml successfully found and loaded')
            tree = ET.parse(lrg_xml)
        except IOError:
            print("Could not find LRG_" + str(sys_args.local) + ".xml locally, it was downloaded instead.")
            url = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_' + str(sys_args.local) + '.xml'
            try:  # Check it worked or throw up an error message
                lrg_xml = urlopen(url)
                tree = ET.parse(lrg_xml)
                tree.write(open('LRG_' + str(sys_args.local) + '.xml', 'wb'))  # Write to file
            except Exception as err:
                print("The file could not be retrieved from the web url - check file name and internet connection?")
                sys.exit(err)  # Exit with an error
    elif sys_args.web:  # Otherwise they want on from the web, so fetch this
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_' + str(sys_args.web) + '.xml'
        try:  # Check it worked or throw up an error message
            lrg_xml = urlopen(url)
            tree = ET.parse(lrg_xml)
            tree.write(open('LRG_' + str(sys_args.web) + '.xml', 'wb'))  # Write to file
        except Exception as err:
            print("The file could not be retrieved from the web url - check file name and internet connection")
            sys.exit(err)  # Exit with an error
    return tree


# cmd_line_input = get_args()
#
# get_lrg_file(cmd_line_input)


def check_version(tree):
    """ Tests that imported XML is an LRG file and has the expected version number """
    root = tree.getroot()
    assert root.tag == 'lrg', 'XML not LRG'
    assert root.get('schema_version') == '1.9', 'Wrong LRG version, should to be 1.9'


def get_summary_list(tree):
    """
    Creates a list and appends it with basic info from the fixed_annotation part of the LRG file.
    Currently picks the lrg_locus, the id (LRG version), the hgnc_id and the sequence source.
    """
    results = list()
    results.append(('Gene', tree.find('updatable_annotation/annotation_set/lrg_locus').text))
    results.append(('LRG version', tree.find('fixed_annotation/id').text))
    results.append(('HGNC ID', tree.find('fixed_annotation/hgnc_id').text))
    results.append(('Sequence source', tree.find('fixed_annotation/sequence_source').text))
    return results


def get_assemblies(tree):
    """
    Finds all LRG assemblies in the annotation_set node and returns a dictionary
    of the assembly name (key) and attributes as nested key/value pair.
    """
    assemblies = {}
    lrg_assembly = tree.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping')
    for item in lrg_assembly:
        assembly_name = item.attrib['coord_system']
        assemblies[assembly_name] = {}
        assemblies[assembly_name]['start'] = item.attrib['other_start']
        assemblies[assembly_name]['end'] = item.attrib['other_end']
        assemblies[assembly_name]['strand'] = item[0].attrib['strand']
    return assemblies


def get_differences(tree):
    """
    Finds all listed differences for each assembly and returns in dictionary of assembly name (key)
    and attributes as nested key/value pair
    """
    differences = {}
    lrg_assembly = tree.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping')
    for assembly in lrg_assembly:
        assembly_name = assembly.attrib['coord_system']
        if assembly.findall('mapping_span/diff'):
            differences[assembly_name] = {}
            diffs = assembly.findall('mapping_span/diff')
            diff_id = 1
            for items in diffs:
                differences[assembly_name][diff_id] = {}
                differences[assembly_name][diff_id]['type'] = items.attrib['type']
                differences[assembly_name][diff_id]['lrg_start'] = items.attrib['lrg_start']
                differences[assembly_name][diff_id]['lrg_end'] = items.attrib['lrg_end']
                differences[assembly_name][diff_id]['lrg_sequence'] = items.attrib['lrg_sequence']
                differences[assembly_name][diff_id]['other_sequence'] = items.attrib['other_sequence']
                diff_id += 1
    return differences


def get_exons(tree):
    """
    Finds all the exons as a list. Returns multiple transcripts if present
    """
    results = {}
    transcripts = tree.findall('fixed_annotation/transcript')
    for ts in transcripts:
        transcript_id = tree.find('fixed_annotation/id').text + ts.attrib['name']
        results[transcript_id] = {}
        for e in ts.findall('exon'):
            for coord in e:
                if (coord.attrib['coord_system']) == transcript_id:
                    results[transcript_id][e.attrib['label']] = (coord.attrib['start'], coord.attrib['end'])
    return results


def output_results(tree):
    """
    Output all results to text file: runs the 'get' functions and iterates over returned data and prints
    """
    results = get_summary_list(tree)

    filename = results[1][1] + ".txt"

    # os.makedirs(os.path.dirname(filename), exist_ok=True)

    os.makedirs('output/', exist_ok=True)

    with open('output/' + filename, 'w') as f:
        f.write("File Summary:\n---------------\n")

        # Output from the get_summary_list function

        for item in results:
            f.write(item[0] + ": " + item[1] + '\n')

        # Output from the get_assemblies function

        f.write("\nLRG Assemblies:\n---------------\n")
        assemblies = get_assemblies(tree)
        f.write("Assembly\t Genomic coordinates\t Strand\n")
        for key1, value1 in assemblies.items():
            f.write(key1 + " \t " + assemblies[key1]['start'] + " - " + assemblies[key1]['end']
                    + " \t " + assemblies[key1]['strand'] + "\n")

        # Output from the get_differences function

        differences = get_differences(tree)
        for key1, value1 in differences.items():
            f.write("\nDifferences with reference Seq:" + key1 + "\n---------------\n")
            f.write("Type\t\tChange\t\tLRG coordinates\n")
            for key2, value2 in value1.items():
                f.write(value2['type'] + " \t" + value2['lrg_sequence'] + ">" + value2['other_sequence'] +
                        "\t\t" + value2['lrg_start'] + " - " + value2['lrg_end'] + "\n")

        # Output from the get_exons function

        f.write("\nExon details:\n---------------\n")
        exons = get_exons(tree)
        for transcript, exons in sorted(exons.items()):
            f.write("Transcript ID: " + transcript + "\n")
            f.write("Exon no.\tLRG coordinates\n")
            for exon_id in sorted(exons.keys(), key=int):
                f.write(exon_id + ":\t\t" + exons[exon_id][0] + "-" + exons[exon_id][1] + "\n")
            f.write("\n")
        print('Results have been saved in output/' + filename)

    f.close()


def main():
    """
    Main function for script - not yet finished.
    """
    sys_args = get_args()
    tree = get_lrg_file(sys_args)
    check_version(tree)
    output_results(tree)


main()  # Call the main function


# if __name__ == '__main__':  # Will sort this out in the final thing
#     main()
