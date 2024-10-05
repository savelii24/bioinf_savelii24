### BioInfUt 

## What is it?
**BioInfUt** is a bioinformatics utility for performing basic RNA/DNA operations and filtering FASTQ sequences.

## Main Features
**BioInfUt** can assist you with various RNA/DNA operations:

- The program accepts one or more DNA or RNA sequences and returns their reversed, complemented, transcribed, or reverse-complemented forms.
- It can also detect sequences that contain both U (Uracil) and T (Thymine) nucleotides in a single molecule.

In addition to these basic RNA/DNA operations, the program can filter FASTQ sequences:

- The program takes a dictionary of FASTQ sequences as input and filters them based on several criteria: GC-content (as a percentage), sequence length (number of nucleotides), and read quality (using the Phred33 scale).

## Where to get it
The source code is hosted on GitHub at:
[https://github.com/savelii24/bioinf_savelii24/bioinf_ut](https://github.com/savelii24/bioinf_savelii24/bioinf_ut)

You can also install the program from the command line:
```sh
# Clone the repository
git clone git@github.com:savelii24/bioinf_savelii24/bioinf_ut.git
``` 
## Also see this

File has the following structure:

   ```python
    - bioinf_ut/
     	|- README.md
     	|- bioinf_ut.py
     	|- modules_4main_script/
           	|- module_for_filter_fastq.py
           	|- module_for_dna_rna_tools.py
           	|- ...
    ``` 
