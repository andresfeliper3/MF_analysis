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

### Analyze and find kmers recursively
Analyze using MFA and find kmers recursively.
This command can be used to execute the whole genome.

    py .\command.py find_kmers_genome -method r -name "caenorhabditis elegans"

This command can be used to execute a single chromosome.

    py .\command.py find_kmers_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -method r -name "caenorhabditis elegans" 

Saving kmers to database is NOT implemented yet.

### Graph and xlsx file
Load and graph the data:

        py .\command.py graph -name "caenorhabditis elegans" -mode whole
        py .\command.py graph -name "caenorhabditis elegans" -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode whole

Compare results with branch main and changeAlgorithm.


## Repeats using RepeatMasker
After executing the RepeatMasker program, the results are saved in a .out file.
These files can be saved in the following directory _resources/RM_resources/<organism_name>_

### Save RepeatMasker results to database 
Save the repeats found by RepeatMasker into the database using the results file path as a parameter.

Example using the command:

    py .\command.py load_RM_repeats -path resources/RM_resources/c_elegans_chromosome_I.fasta.out
