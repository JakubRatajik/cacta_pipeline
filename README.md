# CACTA pipeline

[TOC]

## Introduction

The CACTA pipeline is a tool for structural detection of transposons from the CACTA superfamily at the whole-genome level. CACTA candidates and their annotations are first obtained from the genome using a Python module and then, if required, filtered and sorted into families using a shell script utilizing [VSEARCH](https://github.com/torognes/vsearch)[^1] clustering.

## Prerequisities

- Git for cloning the repository. Alternatively, a [ZIP](https://gitlab.fi.muni.cz/xratajik/cacta_pipeline/-/archive/master/cacta_pipeline-master.zip) file can be downloaded.

- Python 3.7 or newer

- pip - if you have Python installed and added to the $PATH, then simply run

    `python -m ensurepip --upgrade`

- VSEARCH

    - Binaries, source code as well as installation manual are available [here](https://github.com/torognes/vsearch).
    - To clone the repository and build VSEARCH yourself:
      ```
      git clone https://github.com/torognes/vsearch.git
      cd vsearch
      ./autogen.sh
      ./configure CFLAGS="-O3" CXXFLAGS="-O3"
      make
      make install
      ```

</details>`

## Installation

```
# clone the repository
git clone https://gitlab.fi.muni.cz/xratajik/cacta_pipeline
cd cacta_pipeline

# allow to execute the shell script (only relevant to Unix-like OS)
chmod u+x cacta_families.sh

# install required python libraries
pip install -r requirements.txt
```

## Usage

### Detect CACTA

The Detect CACTA module structurally identifies TIRs (Tandem Inverted Repeats) and corresponding TSDs (Target Site duplications) that could potentially dealineate a CACTA transposon. CACTA transposon have characteristic "CACT[A/G]" motif appearing at the distant end of their TIRs.

```
Usage: python detect_cacta.py [options]

Available options:  
  -i [file],            Input DNA sequence file (default: None)
    --in-file [file]     
    
  -fa [file],           Output FASTA file with candidate sequences (default:
    --fasta-out [file]     None)
   
  -g [file],			Output GFF3 file with candidate annotation (default:
    --gff3 [file]           None)
  
  --min-len [int]       Miminum transposon length (default: 50)
  
  --max_len [int]       Maximum transposon length (default: 23018)
  
  --tir-info            Append TIR info (TIR length, mismatch count, gap
                        count) to element name (default: False)
                        
  -h, --help            Prints this help message.
```

### CACTA pipeline

The CACTA pipeline module is a shell-based script employing VSEARCH clustering to filter spurious candidates and assign transposons to families.

```
Usage: ./cacta_families.sh -i <file> [-m <int>] [-e] [-g <file>] [-c] [-h]

Available options:
-i <file>        Input FASTA file with CACTA TE sequences.

-m <int>         Minimum family members. Families with number of
                   members lower than <-m> are discarded.
                   Default value is 2.
                   
-c               FASTA file "<input file>_consensi.fasta" containing
                   consensus sequence for each family is generated.
                   
-e               FASTA file "<input file>_filtered.fasta" containing
                   filtered elements is generated.
                   
-g <file>        GFF3 annotation "<input file>_filtered.gff3" of
                   filtered elements is generated. GFF3 file generated
                   by detect_cacta.py must be provided!
                   
-h               Prints this help message.
```

## Testing

To check if all dependencies were installed correctly, it is recommended to test the pipeline on a test genome first. 

```
python detect_cacta.py -i test/test_genome.fasta -fa test/candidates.fasta -g test/candidates.gff3 --min-len 50 --max-len 23018
./cacta_families.sh -i test/candidates.fasta -g test/candidates.gff3 -m 2 -e -c

```

If everything was set up correctly, the overall analysis should take no more than five minutes. The Detect CACTA should identify 10 candidates and generate a FASTA file containing the element sequences and a GFF3 file containing the annotation.

CACTA families should generate a single family consisting of two members, their annotation and a consensus sequence in respective files. 

---

[^1]: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)