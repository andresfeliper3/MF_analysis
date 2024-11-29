@echo off
setlocal

rem Editable parameters
set NAME="neodiprion lecontei"
set WINDOW_LENGTH=3000000
set DIR="neodiprion_lecontei"
set GCF=GCF_021901455.1
set GTF_PATH=resources/genes/neodiprion_lecontei/gtf/GCF_021901455.1_iyNeoLeco1.1_genomic.gtf
set SEQUENCE_PATH=resources/dna_sequences/neodiprion_lecontei/chr1.fna
set PARTITIONS=300
set REGIONS=3
set K_RANGE="(4,12)"
set GENES_KMER_SIZE=4
set GENES_AMOUNT=100

rem Execute commands
::echo Running download...
::py .\command.py download -name %NAME%

echo Running analyze_sequence -mode whole...
py .\command.py analyze_sequence -path %SEQUENCE_PATH% -name %NAME% -mode whole

echo Running analyze_sequence -mode regions...
py .\command.py analyze_sequence -path %SEQUENCE_PATH% -name %NAME% -mode regions -window_length %WINDOW_LENGTH%

echo Running find_kmers_sequence linearly...
py .\command.py find_kmers_sequence -path %SEQUENCE_PATH% -method l -k_range %K_RANGE% -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR% -graph_from_file true

echo Running graph -mode whole...
py .\command.py graph -name %NAME% -mode whole

echo Running graph -mode region...
py .\command.py graph -name %NAME% -mode region -window_length %WINDOW_LENGTH%

echo Running graph_linear_repeats_sequence...
py .\command.py graph_linear_repeats_sequence --save true -name %NAME% -path %SEQUENCE_PATH% -dir %DIR% -k_range %K_RANGE%

echo Running graph_gtf_file...
py .\command.py graph_gtf_file -path %GTF_PATH% -partitions %PARTITIONS% -regions %REGIONS% -plot_type line -dir %DIR% --save true

echo Running graph_linear_regression_sequence...
py .\command.py graph_linear_regression_sequence -path %SEQUENCE_PATH% -k_range %K_RANGE% -dir %DIR% -name %NAME%

echo Running load_genes...
py .\command.py load_genes -path %GTF_PATH%

:: Last three are dependent
echo Running find_kmers_linearly_genes_sequence...
py .\command.py find_kmers_linearly_genes_sequence -path %SEQUENCE_PATH% -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR% -graph_from_file true -size %GENES_KMER_SIZE%

echo Running graph_linear_in_genes_repeats_sequence
py .\command.py graph_linear_in_genes_repeats_sequence --save true -name %NAME% -path %SEQUENCE_PATH% -dir %DIR% -k_range %K_RANGE%

echo Scraping functional categories and subcategories from KEGG..
py .\command.py load_categories -path %SEQUENCE_PATH% -name %NAME% -size %GENES_KMER_SIZE% -genes_amount %GENES_AMOUNT%

echo Graphing functional categories and subcategories heatmaps...
py .\command.py graph_categories_repeats_heatmap_sequence -path %SEQUENCE_PATH% -size %GENES_KMER_SIZE% --save true -dir %DIR% -name %NAME% --tags true


endlocal