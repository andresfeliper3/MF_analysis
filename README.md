# MF_analysis
# Table of Contents


1. [Configurable Files](#configurable-files)
2. [Commands Examples](#commands-examples)
   - [Download FASTA Sequences](#download-fasta-sequences)
   - [Analyze](#analyze)
   - [Analyze and Find kmers Recursively](#analyze-and-find-kmers-recursively)
3. [Executing RepeatMasker](#executing-repeatmasker)
   - [RepeatMasker for a Single Sequence](#repeatmasker-for-a-single-sequence)
   - [RepeatMasker Execution for Genome Folder](#repeatmasker-execution-for-genome-folder)
4. [Graph](#graph)
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
5. [Repeats Using RepeatMasker](#repeats-using-repeatmasker)
   - [Save RepeatMasker Results to Database](#save-repeatmasker-results-to-database)
     - [Using a Folder with the Results of a Genome](#using-a-folder-with-the-results-of-a-genome)
6. [Genes Using a .gtf File](#genes-using-a-gtf-file)
   - [Save the .gtf File Data to Database](#save-the-gtf-file-data-to-database)


## Configurable files
* _resources/sequences.yaml_ contains the chromosomes info to execute the MFA analysis. It can be modified following the predefined format.
The organism name and GCF expressed in that file are the ones to be used in the commands.

* _resources/db_config.yaml_ contains the configuration of the selected database to work with. Modify the *type* value in order to select the type of database.

## Commands examples

### Download FASTA sequences
Download and uncompress files

        py .\command.py download -name "caenorhabditis elegans"
        py .\command.py download -name GCF_000002985.6

### Analyze
Analyze and load the organism and genome data. 

                    
        py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions
        py .\command.py analyze_genome -name GCF_000002985.6 -mode regions -window_length 300000
        py .\command.py analyze_genome -name GCF_000002985.6 -mode regions -regions_number 3
        py .\command.py analyze_genome -name GCF_000002985.6 -mode whole

Analyze only one sequence file (one chromosome) given a file path.
When analyzing by regions, choose between dividing the sequence by amount of regions or by region length.
Use one of these parameters:
- window_length: size of the region/window in bps. The amount of windows will be deduced.
- regions_number: amount of regions that will be used (all of them of the same size except for the last one). The window length is calculated.
    

        py .\command.py analyze_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name GCF_000002985.6 -mode regions -window_length 300000
        py .\command.py analyze_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name GCF_000002985.6 -mode regions -regions_number 3 
        py .\command.py analyze_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name GCF_000002985.6 -mode whole

It is possible to configure if the results should be saved to the database or not.
The default behavior is to always save the data to database.

For the whole genome:

         py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions --save-to-db=false
         py .\command.py analyze_genome -name "caenorhabditis elegans" -mode regions --save-to-db=true

For a single sequence:

         py .\command.py analyze_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -mode whole --save-to-db=false
         py .\command.py analyze_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -mode whole --save-to-db=true

### Analyze and find kmers recursively
Analyze using MFA and find kmers recursively.
This command can be used to execute the whole genome.

    py .\command.py find_kmers_genome -method r -k_range "(4, 8)" -name "caenorhabditis elegans"

This command can be used to execute a single chromosome.

    py .\command.py find_kmers_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -method r  -k_range (4, 8) -name "caenorhabditis elegans" 


### Analyze and find kmers linearly
Find kmers going through the sequence linearly.
If graph_from_file is true, then the graphs are generated from reading the files instead of using the database.
In this case, adding window_length a value is mandatory for a good performance.

- -method: l
- -k_range: "(4, 12)"
- -window_length: 300000
- -graph_from_file: true

This command can be used to execute the whole genome.

    py .\command.py find_kmers_genome -method l -k_range "(4, 12)" -name "caenorhabditis elegans" -window_length 300000 -dir "caenorhabditis_elegans" -graph_from_file true

This command can be used to execute a single chromosome
    
    py .\command.py find_kmers_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -method l -k_range "(4, 8)" -name "caenorhabditis elegans" -window_length 300000 -dir "caenorhabditis_elegans" -graph_from_file true

### Analyze and find kmers linearly only for kmers found in genes
Find the kmers linearly. Find the most common kmers and get their frequency only in genes.

- -path: sequence path file.
- -name: organism name in sequence yaml file.
- -window_length
- -dir: subfolder name to store the graphs.
- -graph_from_file: (true) generate graphs directly from sequence files.
- -size: size of the kmers/repeats that are going to be searched in the genes.

    
    py .\command.py find_kmers_linearly_genes_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -window_length 300000 -dir "caenorhabditis_elegans" -graph_from_file true -size 4

Execution for genomes:
- -name: it helps identify the folder where all the chromosomes files are.


    py .\command.py find_kmers_linearly_genes_genome -name "caenorhabditis elegans" -window_length 300000 -dir "caenorhabditis_elegans" -graph_from_file true -size 4

## Executing RepeatMasker
In order to use RepeatMasker from a Docker container. 
Execute the container using the following command:

    docker run -it --rm --mount type=bind,source="D:/TG/MFA_analysis/resources",target=/working --mount type=bind,source="E:/lib",target=/opt/RepeatMasker/Libraries/famdb dfam/tetools:latest bash 

The first *source* is the path to the /resources folder of the project. 
The second *source* is the path to the /lib folder where all the partitions of the Dfam database are located.

### RepeatMasker for a single sequence
In the Docker container, go to the /working directory and execute the following command:

    RepeatMasker -species nematode -dir /working/RM_resources/caenorhabditis_elegans /working/dna_sequences/caenorhabditis_elegans/c_elegans_chromosome_X.fasta

You can change the *-dir* parameter to change the destination folder. The next parameter is the path to the .fasta or 
.fna sequence.

### RepeatMasker execution for genome folder
In the Docker container, go to the /working directory and execute the following command. 

    ./repeatmasker_species_runner.sh /working/dna_sequences/Acropora_millepora/ /working/RM_resources/acropora_millepora/ cnidaria

The first parameter is the *input directory* and the second parameter is the *output directory*.
Change the input directory organism name <Acropora_millepora> in order to change the organism to be executed.
Change the output directory organism name <acropora_millepora> in order to change the output folder of the results.
Remember to have all the chromosomes files .fna or .fasta in the input directoryy.

## Graph

### Graph and xlsx file
Load and graph the data:

        py .\command.py graph -name "caenorhabditis elegans" -mode whole
        py .\command.py graph -name GCF_000002985.6 -mode regions -regions_number 3 
        py .\command.py graph -name GCF_000002985.6 -mode regions -window_length 300000
        py .\command.py graph -name GCF_000002985.6 -mode whole

Compare results with branch main and changeAlgorithm.

### Graph the RepeatMasker results from a RM results file
The graphs are saved in the /out directory in the sequence folder.

It is mandatory to specify the
- path of the RM results file.
- dir - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
- partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_rm_file_sequence -path resources/RM_resources/caenorhabditis_elegans/c_elegans_chromosome_I.fasta.out -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_file_sequence -path resources/RM_resources/caenorhabditis_elegans/c_elegans_chromosome_I.fasta.out -dir caenorhabditis_elegans
 
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

         py .\command.py graph_rm_file_genome -path resources/RM_resources/caenorhabditis_elegans -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_file_genome -path resources/RM_resources/caenorhabditis_elegans -dir caenorhabditis_elegans 


### Graph the Plantrep.cn repeats using the file
The graphs are saved in the /out directory in the sequence folder. 
Keep in mind that the repeats file contains the repeats of the entire genome of the organism,
not only a single chromosome. 

It is mandatory to specify the
- path of the repeats file.
- dir - scientific name of the organism that will be used as the folder name to save the graph.

Other possible parameters are:
  - partitions - number of partitions to use in the graph generation. Each partition will represent one point in the graph.
    300 by default.
  - regions - amount of regions in which the graph will be divided using vertical line. 3 by default
  - plot_type - style of plot (line or bar), line is by default.
  - save - save the graph in the local directory specified /out (true or false).

         py .\command.py graph_genome_repeats_from_file -path resources/RM_resources/musa_acuminata/157_Musa_acuminata_rm.out -partitions 300 -regions 3 -plot_type line -dir musa_acuminata --save true
        
The command without the optional parameters would be:

       py .\command.py graph_genome_repeats_from_file -path resources/RM_resources/musa_acuminata/157_Musa_acuminata_rm.out -dir musa_acuminata
 


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

         py .\command.py graph_rm_database_sequence -ran NC_003279.8 -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true

The command without the optional parameters would be:

       py .\command.py graph_rm_database_sequence -ran NC_003279.8 -dir caenorhabditis_elegans

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

         py .\command.py graph_rm_database_genome -gcf GCF_000002985.6 -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true
        
The command without the optional parameters would be:

       py .\command.py graph_rm_database_genome -gcf GCF_000002985.6 -dir caenorhabditis_elegans


### Graph the recursively found repeats results from the database
The graphs are saved in the /out directory in the sequence folder.

#### Graph per sequence/chromosome
The n_max parameter is optional. It represents the total amount of repeats shown in the general graph.
The chromosome is identified by the refseq_accession_number (-ran)

    py .\command.py graph_recursive -ran NC_003279.8 --save true -name "caenorhabditis_elegans" 
    py .\command.py graph_recursive -ran NC_003279.8 --save true -name "caenorhabditis_elegans" -n_max 10

#### Graph per genome
The n_max parameter is optional. It represents the total amount of repeats shown in the general graphs.
The genome is identified by the GCF (-gcf).

    py .\command.py graph_recursive_genome -gcf GCF_000002985.6 --save true -dir caenorhabditis_elegans 
    py .\command.py graph_recursive_genome -gcf GCF_000002985.6 --save true -dir caenorhabditis_elegans -n_max 10

### Graph repeats found linearly from the database
Graph the repeats found with the linear algorithm using a sequence path and results saved in the database.

- --save: true (to save graphs)
- -name: name of the organism to read it from yaml file.
- -path: relative path to the sequence.
- -window_length: window size for region analysis.
- -dir: directory to save the graphs (organism_name).
- -k_range: range of values of k in kmers.


    py .\command.py graph_linear_repeats_sequence --save true -name "Caenorhabditis elegans" -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -dir caenorhabditis_elegans -k_range "(4,12)"


To graph the entire genome use the following command:

    py .\command.py graph_linear_repeats_genome --save true -name "Caenorhabditis elegans" -gcf GCF_000002985.6 -dir caenorhabditis_elegans -k_range "(4,12)

### Graph repeats found linearly only in genes from the database
Graph the repeats found with the linear algorithm ONLY IN GENES using a sequence path and results saved in the database.

- --save: true (to save graphs)
- -name: name of the organism to read it from yaml file.
- -path: relative path to the sequence.
- -window_length: window size for region analysis.
- -dir: directory to save the graphs (organism_name).
- k_range: the range of k values to consider when graphing. Ex: "(4, 12)".


    py .\command.py graph_linear_in_genes_repeats_sequence --save true -name "Caenorhabditis elegans" -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -dir caenorhabditis_elegans -k_range "(4,12)"

If you want to run the entire genome:
Instead of path, use gcf:
- -gcf: organism GCF.


    py .\command.py graph_linear_in_genes_repeats_genome --save true -name "Caenorhabditis elegans" -gcf GCF_000002985.6 -dir caenorhabditis_elegans -k_range "(4,12)"


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

    py .\command.py graph_gtf_file -path  resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true

Short version:

    py .\command.py graph_gtf_file -path  resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf -dir caenorhabditis_elegans

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

    py .\command.py graph_gtf_database -gcf GCF_000002985.6 -partitions 300 -regions 3 -plot_type line -dir caenorhabditis_elegans --save true

Short version:

    py .\command.py graph_gtf_database -gcf GCF_000002985.6 -dir caenorhabditis_elegans
    py .\command.py graph_gtf_database -ran NC_003279.8  -dir caenorhabditis_elegans 
    py .\command.py graph_gtf_database -ran NC_003279.8  -dir caenorhabditis_elegans -gcf anything_here

## Graph linear regression
### Linear regression DDq vs repeats frequency in regions
- -path: path of the sequence.
- -k_range
- -dir: directory to save the graphs.
- name: name of the organism name.


    py .\command.py graph_linear_regression_sequence -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -k_range "(4,12)" -dir caenorhabditis_elegans -name "caenorhabditis elegans"

For the genome:

    py .\command.py graph_linear_regression_genome -gcf GCF_000002985.6 -k_range "(4,12)" -dir caenorhabditis_elegans -name "caenorhabditis elegans"


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

#### Using a file with the results of a genome
The file is segmented in chromosomes before saving it to the database.
Example using the command:

    py .\command.py load_genome_repeats_file -path resources/RM_resources/musa_acuminata/157_Musa_acuminata_rm.out


## Genes using a .gtf file
### Save the .gtf file data to database 
Example using the command:

    py .\command.py load_genes -path resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf 

## Add categories and subcategories from KEGG website
Use this command to scrape the KEGG website and link the identified genes for the repeats 
with their categories and subcategories of KEGG pathway.
- -path: path of the sequence that whose genes are going to be linked to KEGG categories and subcategories.
- -name: name of the organism in YAML file.
- -size: size of the kmers/repeats whose genes are going to be linked with KEGG info.


    py .\command.py load_categories -path resources/dna_sequences/caenorhabditis_elegans/chrI.fna -name "caenorhabditis elegans" -size 4

For genome:

    py .\command.py load_categories_genome -name "caenorhabditis elegans" -size 4



