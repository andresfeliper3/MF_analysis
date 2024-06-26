# MF_analysis


## Configurable files
The file _resources/sequences.yaml_ contains the chromosomes info to execute the MFA analysis. It can be modified following the predefined format.
The organism name and GCF expressed in that file are the ones to be used in the commands.

## Commands examples

### Download FASTA sequences
Download and uncompress files

        py .\command.py download -name "caenorhabditis elegans"
        py .\command.py download -name GCF_000002985_4

### Analyze
Analyze and load the organism and genome data. 

                    
        py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions
        py .\command.py analyze_genome -name GCF_000002985_4 -mode regions
        py .\command.py analyze_genome -name GCF_000002985_4 -mode whole

Analyze only one sequence file (one chromosome) given a file path
    
        py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -mode regions
        py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name GCF_000002985_4 -mode regions
        py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name GCF_000002985_4 -mode whole

It is possible to configure if the results should be saved to the database or not.
The default behavior is to always save the data to database.

For the whole genome:

         py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions --save-to-db=false
         py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions --save-to-db=true

For a single sequence:

         py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -mode whole --save-to-db=false
         py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -mode whole --save-to-db=true


### Graph and xlsx file
Load and graph the data:

        py .\command.py graph -name "caenorhabditis elegans" -mode whole
        py .\command.py graph -name "caenorhabditis elegans" -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode whole

Compare results with branch main and changeAlgorithm