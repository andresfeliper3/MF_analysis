# MF_analysis
# Table of Contents

1. [Configurable Files](#configurable-files)
2. [Commands Examples](#commands-examples)
   - [Download FASTA Sequences](#download-fasta-sequences)
   - [Analyze](#analyze)
   - [Analyze and Find kmers Recursively](#analyze-and-find-kmers-recursively)
3. [Graph](#graph)
   - [Graph and XLSX File](#graph-and-xlsx-file)
   - [Graph the RepeatMasker Results from a RM Results File](#graph-the-repeatmasker-results-from-a-rm-results-file)
   - [Graph Using a Genome Folder of .out Result Files](#graph-using-a-genome-folder-of-out-result-files)
   - [Graph Using a RefSeq Accession Number and the Database](#graph-using-a-refseq-accession-number-and-the-database)
   - [Graph the Recursively Found Repeats Results from the Database](#graph-the-recursively-found-repeats-results-from-the-database)
     - [Graph per Sequence/Chromosome](#graph-per-sequencechromosome)
     - [Graph per Genome](#graph-per-genome)
   - [Graph Genes Data](#graph-genes-data)
     - [Graph from .gtf File](#graph-from-gtf-file)
     - [Graph from Database](#graph-from-database)
4. [Repeats using RepeatMasker](#repeats-using-repeatmasker)
   - [Save RepeatMasker Results to Database](#save-repeatmasker-results-to-database)
     - [Using a Folder with the Results of a Genome](#using-a-folder-with-the-results-of-a-genome)
5. [Genes using a .gtf File](#genes-using-a-gtf-file)
   - [Save the .gtf File Data to Database](#save-the-gtf-file-data-to-database)


## Configurable files
The file _resources/sequences.yaml_ contains the chromosomes info to execute the MFA analysis. It can be modified following the predefined format.
The organism name and GCF expressed in that file are the ones to be used in the commands.

## Commands examples

### Download FASTA sequences
Download and uncompress files

        py .\command.py download -name "caenorhabditis elegans"
        py .\command.py download -name GCF_000002985_6

### Analyze
Analyze and load the organism and genome data. 

                    
        py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions
        py .\command.py analyze_genome -name GCF_000002985_6 -mode regions
        py .\command.py analyze_genome -name GCF_000002985_6 -mode whole

Analyze only one sequence file (one chromosome) given a file path
    

        py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name GCF_000002985_6 -mode regions
        py .\command.py analyze_sequence -path resources/dna_sequences/Caenorhabditis_elegans/chrI.fna -name GCF_000002985_6 -mode whole

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

    py .\command.py find_kmers_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -method r -name "caenorhabditis elegans" 

Saving kmers to database is NOT implemented yet.

## Graph

### Graph and xlsx file
Load and graph the data:

        py .\command.py graph -name "caenorhabditis elegans" -mode whole
        py .\command.py graph -name "caenorhabditis elegans" -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode regions
        py .\command.py graph -name GCF_000002985_4 -mode whole

Compare results with branch main and changeAlgorithm.

### Graph the RepeatMasker results from a RM results file
The graphs are saved in the /out directory in the sequence folder.

It is mandatory to specify the
- path of the RM results file.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_rm_file_sequence -path resources/RM_resources/caenorhabditis_elegans/c_elegans_chromosome_I.fasta.out -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_file_sequence -path resources/RM_resources/caenorhabditis_elegans/c_elegans_chromosome_I.fasta.out -name "Caenorhabditis elegans" 
 
#### Graphing using a genome folder of .out result files 
It is mandatory to specify the
- path of the RM results folder.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_rm_file_genome -path resources/RM_resources/caenorhabditis_elegans -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_file_genome -path resources/RM_resources/caenorhabditis_elegans -name "Caenorhabditis elegans" 


### Graph the RepeatMasker results using a refseq accession number and the database
It is mandatory to specify the
- ran - refseq accession number, to identify the sequence.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_rm_database_sequence -ran NC_003279.8 -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true

The command without the optional parameters would be:

       py .\command.py graph_rm_database_sequence -ran NC_003279.8 -name "Caenorhabditis elegans" 

#### Graphing using a genome GCF and the database
It is mandatory to specify the
- gcf - GCF to identity the organism genome.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_rm_database_genome -gcf GCF_000002985.6 -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_database_genome -gcf GCF_000002985.6 -name "Caenorhabditis elegans" 


### Graph the recursively found repeats results from the database
The graphs are saved in the /out directory in the sequence folder.

#### Graph per sequence/chromosome
The n_max parameter is optional. It represents the total amount of repeats shown in the general graph.
The chromosome is identified by the refseq_accession_number (-ran)

    py .\command.py graph_recursive -ran NC_003279.8 --save true -name "caenorhabditis elegans" 
    py .\command.py graph_recursive -ran NC_003279.8 --save true -name "caenorhabditis elegans" -n_max 10

#### Graph per genome
The n_max parameter is optional. It represents the total amount of repeats shown in the general graphs.
The genome is identified by the GCF (-gcf).

    py .\command.py graph_recursive_genome -gcf GCF_000002985.6 --save true -name "caenorhabditis elegans" 
    py .\command.py graph_recursive_genome -gcf GCF_000002985.6 --save true -name "caenorhabditis elegans" -n_max 10

### Graph genes data
#### Graph from .gtf file
It is mandatory to specify the
- path - relative path of the .gtf file.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

Example using the command:

    py .\command.py graph_gtf_file -path  resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true

Short version:

    py .\command.py graph_gtf_file -path  resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf -name "Caenorhabditis elegans" 

#### Graph from database
It is mandatory to specify the
- gcf - GCF of the organism/genome.
  - ran - refseq accession number of the sequence/chromosome.
  - name - scientific name of the organism that will be used as the folder name to save the graph.

**If GCF is entered, it is NOT necessary to enter a RAN. If RAN is entered, it is NOT necessary to enter
a GCF. If RAN is entered, any GCF added will be ignored.**

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

Example using the command:

    py .\command.py graph_gtf_database -gcf GCF_000002985.6 -partitions 300 -regions 3 -plot_type line -name "Caenorhabditis elegans" --save true

Short version:

    py .\command.py graph_gtf_database -gcf GCF_000002985.6 -name "Caenorhabditis elegans"  
    py .\command.py graph_gtf_database -ran NC_003279.8  -name "Caenorhabditis elegans"
    py .\command.py graph_gtf_database -ran NC_003279.8  -name "Caenorhabditis elegans" -gcf anything_here


## Repeats using RepeatMasker
After executing the RepeatMasker program, the results are saved in a .out file.
These files can be saved in the following directory _resources/RM_resources/<organism_name>_

### Save RepeatMasker results to database 
Save the repeats found by RepeatMasker into the database using the results file path as a parameter.

Example using the command:

    py .\command.py load_RM_repeats -path resources/RM_resources/caenorhabditis_elegans/c_elegans_chromosome_I.fasta.out

#### Using a folder with the results of a genome
Example using the command:

    py .\command.py load_RM_repeats_folder -path resources/RM_resources/caenorhabditis_elegans


## Genes using a .gtf file
### Save the .gtf file data to database 
Example using the command:

    py .\command.py load_genes -path resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf 

