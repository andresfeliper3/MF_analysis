@echo off
setlocal

rem Editable parameters
set NAME="saccharomyces mikatae"
set WINDOW_LENGTH=30000
set DIR="saccharomyces_mikatae"
set GCF=GCF_947241705.1
set GTF_PATH=resources/genes/saccharomyces_mikatae/gtf/GCF_947241705.1_Smik-IFO1815_genomic.gtf
set PARTITIONS=100
set REGIONS=3
set K_RANGE="(4,12)"
set GENES_KMER_SIZE=4
set GENES_AMOUNT=100

rem Execute commands

echo Running graph -mode whole...
py .\command.py graph -name %NAME% -mode whole

echo Running graph -mode regions...
py .\command.py graph -name %NAME% -mode regions -window_length %WINDOW_LENGTH%

echo Running graph_linear_repeats_genome...
py .\command.py graph_linear_repeats_genome --save true -name %NAME% -gcf %GCF% -dir %DIR% -k_range %K_RANGE%

echo Running graph_gtf_file...
py .\command.py graph_gtf_file -path %GTF_PATH% -partitions %PARTITIONS% -regions %REGIONS% -plot_type line -dir %DIR% --save true

echo Running graph_linear_regression_genome...
py .\command.py graph_linear_regression_genome -gcf %GCF% -k_range %K_RANGE% -dir %DIR% -name %NAME%

echo Running Graphing functional categories and subcategories heatmaps...
py .\command.py graph_categories_repeats_heatmap_genome -size %GENES_KMER_SIZE% --save true -dir %DIR% -name %NAME% --tags true


endlocal