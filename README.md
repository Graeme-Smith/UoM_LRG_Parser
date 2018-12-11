#lrgParser

---

#### A LRG parser by Graeme Smith and Callum Rakhit, created for the BIOL68400 Bioinformatics Programming Assignment

`lrgParser.py` is a [Python3](https://www.python.org/download/releases/3.0/) script which parses LRG (XML) files and saves the extracted information as a .txt and .BED file. 

The .txt file contains details on the gene (HGNC naming format), LRG version, HGNC ID, sequence source, LRG assembly version, exons, and differences between the genomic reference sequences.

The .BED file contains the chromosome number, genomic coordinates (GRCh38.p12) and the strand of a given sequence.

It can be used to parse a local LRG file:

`python3 Lab_Activity_Data.R -l path/to/your/local/lrg_file.xml`

Also note that if the LRG_ID.xml file you specified is not found the relevant LRG_ID.xml will be downloaded instead, e.g. LRG_1.xml.

The tool can download the relevant LRG file based on a HGNC gene name or LRG file ID:

`python3 Lab_Activity_Data.R -g EGFR`

Created on 21st November 2018.

Please see the documentation.docx for more a more detailed explanation on the usage of this tool.
