# Bioinformatics BIOL68400 Programming Assignment: Parser for LRG (XML) files
# Authors: Callum Rakhit and Graeme Smith
# Date: 21 November 2018
# Version 0.1
#
# This file: LRG_parser.py:
# Script downloads a XML from the LRG website, parses it with ElementTree and prints out gene information
# such as LRG start and end base positions of each exon which is then be used to create the BED file (currently
# just a .txt file).

# Import libraries

from urllib.request import urlopen  # Python 3 specific
import xml.etree.ElementTree as ElTr
import sys
import argparse
import os


def get_args():
    """
    Function uses argparse to setup command line arguments and catch errors
    """
    parser = argparse.ArgumentParser()
    file_location = parser.add_mutually_exclusive_group(required=True)
    file_location.add_argument('-l', '--local', type=int, action='store',
                               help='Takes the LRG ID and parses a copy of the lrg'
                                      'file from the local directory. Assumes file is using '
                                      'the same naming convention as the LRG website. '
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
            tree = ElTr.parse(lrg_xml)
        except IOError:
            print("Cannot find LRG_" + str(sys_args.local) + ".xml, use --web to download")
    elif sys_args.web:  # Otherwise they want on from the web, so fetch this
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_' + str(sys_args.web) + '.xml'
        try:  # Check it worked or throw up an error message
            lrg_xml = urlopen(url)
            tree = ElTr.parse(lrg_xml)
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


#####################################################

