"""
Bioinformatics BIOL68400 Programming Assignment: Parser for LRG (XML) files
Authors: Callum Rakhit and Graeme Smith
Created: 21st November 2018
Description: Parses LRG XML files and outputs a BED file
Usage: See README and documentation.docx for detailed documentation. Made for Python 3.6.6
"""

# TODO
#  - Tidy up function that pulls HGNC positions from the LRG file
#  - Expand .BED file to include strand
#  - Ensure there are adequate unit tests
#  - Ensure there is adequate documentation

import argparse
import os
import requests
import sys
import xml.etree.ElementTree as ElTr

# Import Arguments from command line


def parser_args():
    parser = argparse.ArgumentParser(
        description='Downloads and parses Locus Reference Genomic (LRG) files and produces a .txt and .BED file')
    file_location = parser.add_mutually_exclusive_group(required=False)
    file_location.add_argument(
        '-l', '--local', action='store_true',
        help='Takes a copy of the lrg file from the local directory instead of from the LRG FTP site. Assumes file is '
             'using the same naming convention as the LRG website,  i.e. LRG_{user input}.xml. If the file cannot be '
             'found it is downloaded from the internet instead',
        required=False)
    parser.add_argument(
        '-hgnc', '-g', nargs='+',
        help='Import LRG files for conversion into a BED file as per provided HGNC IDs',
        required=False)
    parser.add_argument(
        '-xref', '-x', nargs='+',
        help='Import LRG files for conversion into a BED file as per provided external references',
        required=False)
    parser.add_argument(
        '-output_dir', '-o', type=str,
        help='File path to output BED file. Defaults to current working directory',
        required=False)
    parser.add_argument(
        '-bed_file', '-b', type=str,
        help='Specify the name of the BED file which the script will output. If not specified will automatically '
             'append .BED file suffix',
        required=False)
    return parser.parse_args()


# Get the LRG ID for an inputted HGNC gene name or external cross reference (Ensembl)


def get_lrg_id(ref, lrg_id_type):
    """
    Uses the EMBL-EBI RESTful API to translate HGNC symbols or external refs into a LRG file name
    """

    # Parse the input arguments, search the web, get LRG file name

    lrg_id = 0

    if lrg_id_type == "hgnc":
        url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=name:' + ref
        response = requests.get(url, allow_redirects=True)
        root = ElTr.fromstring(response.content)
        for entry in root.iter('entry'):
            lrg_id = entry.attrib["id"]

    elif lrg_id_type == "xref":
        url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=' + ref
        response = requests.get(url, allow_redirects=True)
        root = ElTr.fromstring(response.content)
        for entry in root.iter('entry'):
            lrg_id = entry.attrib["id"]

    return lrg_id


# Find the correct LRG file an inputted HGNC gene name or requested LRG number (e.g. LRG_1)


def get_lrg_file(lrg_id, local=False):
    """
    This pulls the LRG/XML file either locally or from the LRG FTP site and
    parses using ElementTree. The "local" flag means a local file is present, if false then
    retrieve file from the LRG website.
    """
    url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/' + lrg_id + '.xml'  # requests doesn't support ftp?

    if local == 1:  # If the user specified a local file
        try:  # Check it worked or throw up an error message
            lrg_xml = requests.get(url, allow_redirects=True)  # Does not support XML, so need a few fixes
            tree = ElTr.fromstring(lrg_xml.content)  # requests function returning a 200 response, changed to str input
            # Now an element not an object, needs to be an object to write to file
            tree = ElTr.ElementTree(tree)  # now an object
            tree.write(open(lrg_id + '.xml', 'wb'))  # write to file
        except IOError:
            print("Could not find " + lrg_id + ".xml locally, it was downloaded instead.")
            try:  # Check it worked or throw up an error message
                lrg_xml = requests.get(url, allow_redirects=True)
                tree = ElTr.fromstring(lrg_xml.content)  # requests was returning a response so changed to string input
                # Now an element not an object, needs to be an object to write to file
                tree = ElTr.ElementTree(tree)  # now an object
                tree.write(open(lrg_id + '.xml', 'wb'))  # Write to file
            except Exception as err:
                print("The file could not be retrieved from the web url - check file name and internet connection?")
                sys.exit(err)  # Exit with an error
            return tree

    elif local == 0:  # Otherwise they want on from the web, so fetch this
        try:  # Check it worked or throw up an error message
            lrg_xml = requests.get(url, allow_redirects=True)
            tree = ElTr.fromstring(lrg_xml.content)  # requests was returning a response so changed to string input
            # Now an element not an object, needs to be an object to write to file
            # tree = ElTr.parse(lrg_xml)  # could previously just parse
            tree = ElTr.ElementTree(tree)  # now an object
            tree.write(open(lrg_id + '.xml', 'wb'))  # Write to file
        except Exception as err:
            print("The file could not be retrieved from the web url - check file name and internet connection")
            sys.exit(err)  # Exit with an error
        return tree


# Check the LRG version is correct


def check_version(tree):
    """ Tests that imported XML is an LRG file and has the expected version number """
    root = tree.getroot()
    assert root.tag == 'lrg', 'XML not LRG'
    assert root.get('schema_version') == '1.9', 'Wrong LRG version, should to be 1.9'


# Functions to populate the .txt or .BED output


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


def get_chromosome(tree):
    """
    Finds chromosome number, start and end positions in the LRG file. Returns them as dictionary of dictionaries.
    """
    chrom_dict = {}
    root = tree.getroot()

    # Creates a dictionary for each set of chromosome co-ordinates, collated in a dictionary

    for mapping in root.iter('mapping'):
        mapping_attributes = mapping.attrib
        chrom_dict[mapping_attributes.get("coord_system", "none")] = mapping_attributes

    # Get the strand for the GRCh38.p12 coord_system

    element = root.find("updatable_annotation/annotation_set/mapping[@coord_system='GRCh38.p12']/mapping_span")
    strand_id = element.get("strand")

    # Hard coded to find positions only from the GRCh38.p12 assembly

    chromosome = list()
    chromosome.append(('chr' + chrom_dict['GRCh38.p12']['other_name']))
    chromosome.append((chrom_dict['GRCh38.p12']['other_start']))
    chromosome.append((chrom_dict['GRCh38.p12']['other_end']))
    chromosome.append(strand_id)  # Looks for the strand

    return chromosome


# Write desired functions to a .txt output


def output_results(tree):
    """
    Output all results to .txt file: runs the 'get' functions and iterates over returned data and prints
    """
    results = get_summary_list(tree)
    filename = results[1][1] + ".txt"
    os.makedirs('output/', exist_ok=True)  # Make the output directory if it doesn't exist

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


# Write desired functions to a .BED output


def output_bed(tree):
    """
    Output all results to a BED file: runs the 'get' functions and iterates over returned data and prints to file
    """
    results = get_summary_list(tree)
    filename = results[1][1] + ".BED"
    chromosome = get_chromosome(tree)
    os.makedirs('output/', exist_ok=True)  # Make the output directory if it doesn't exist

    with open('output/' + filename, 'w') as f:

        # Output from the get_chromosome function

        for item in chromosome:
            f.write(item + '\t')

    f.close()


def main():

    # Parse commandline args

    args = parser_args()

    # Check if HGNC symbols have been specified by user:

    tree = 0

    if args.hgnc is not None:
        for hgnc_ref in args.hgnc:
            tree = (get_lrg_file(get_lrg_id(hgnc_ref, "hgnc"), local=False))

    # Check if external refs have been specified by user:

    elif args.xref is not None:
        for xref_ref in args.xref:
            tree = get_lrg_file(get_lrg_id(xref_ref, "xref"), local=False)

    # Check the version of the LRG file and output function results

    check_version(tree)
    output_results(tree)
    output_bed(tree)


if __name__ == '__main__':  # if this .py script is executing as the main function, the run main()
    main()
