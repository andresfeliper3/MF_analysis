@echo off
setlocal

rem Editable parameters
set NAME="fusarium oxysporum"
set WINDOW_LENGTH=100000
set DIR="fusarium_oxysporum"
set GCF=GCF_000149955.1
set GTF_PATH=resources/genes/fusarium_oxysporum/gtf/GCF_000149955.1_ASM14995v2_genomic.gtf
set PARTITIONS=300
set REGIONS=3
set K_RANGE="(4,12)"
set GENES_KMER_SIZE=4
set GENES_AMOUNT=100

rem Execute commands
::echo Running download...
::py .\command.py download -name %NAME%

echo Running analyze_genome -mode whole...
::py .\command.py analyze_genome -name %NAME% -mode whole

echo Running analyze_genome -mode regions...
py .\command.py analyze_genome -name %NAME% -mode regions -window_length %WINDOW_LENGTH%

echo Running graph -mode whole...
py .\command.py graph -name %NAME% -mode whole

echo Running graph -mode regions...
py .\command.py graph -name %NAME% -mode regions -window_length %WINDOW_LENGTH%

echo Running find_kmers_genome linearly...
py .\command.py find_kmers_genome -method l -k_range %K_RANGE% -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR%

echo Running graph_linear_repeats_genome...
py .\command.py graph_linear_repeats_genome --save true -name %NAME% -gcf %GCF% -dir %DIR% -k_range %K_RANGE%

::echo Running graph_gtf_file...
::py .\command.py graph_gtf_file -path %GTF_PATH% -partitions %PARTITIONS% -regions %REGIONS% -plot_type line -dir %DIR% --save true

::echo Running graph_linear_regression_genome...
::py .\command.py graph_linear_regression_genome -gcf %GCF% -k_range %K_RANGE% -dir %DIR% -name %NAME%

echo Running load_genes...
py .\command.py load_genes -path %GTF_PATH%

:: Last three are dependent
echo Running find_kmers_linearly_genes_genome...
py .\command.py find_kmers_linearly_genes_genome -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR% -graph_from_file true -size %GENES_KMER_SIZE%

echo Running scraping functional categories and subcategories from KEGG..
py .\command.py load_categories_genome -name %NAME% -size %GENES_KMER_SIZE% -genes_amount %GENES_AMOUNT%

echo Running Graphing functional categories and subcategories heatmaps...
py .\command.py graph_categories_repeats_heatmap_genome -size %GENES_KMER_SIZE% --save true -dir %DIR% -name %NAME% --tags true


endlocal