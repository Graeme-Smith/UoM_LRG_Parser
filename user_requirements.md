# LRG Parser
## Participants
- Product owner: Graeme Smith and Callum Rakhit
- Stakeholders: Andy Brass, University of Manchester

## Current Status
Draft
- 21/11/2018

## Purpose
To efficiently pull currated information from an LRG (Locus Reference Genomic) file and convert into a BED file which can seamlessly be used in a wide-range of genomic pipelines.

## Project Goals & Objectives
The tool must satisfied good scientific programming best standards and conform with relevant UKAS, ISO standards, and in-house programming guides.

## Requirements
### Functional
Tool will automate the creation/curation of BED files improving efficiency, removing the risk of manual errors, and allowing faster uptake of information of new transcript information into genomic pipelines.
- Tool will be used by a CLinical Bioinformatician as per a specification for a BED file from a Clinical Scientist.
- A log file detailing creation of the BED file should be produced so that the full providence of a BED file can be retrieved at a later date for troubleshooting purposes.

### Technical (non-functional)
To conform with ISO, UKAS, and in-house standards:
- Code must be well commented and make use of docstrings.
- Any Python code should follow PEP8 style format.
- A README.md file should be present and provide guidance on using the tool.
- Critical code should be covered by unit tests.
- Tool should produce log files detailing the steps taken to produce the BED file.

### Usability
- Tool will have a commandline interface.
