@echo off
setlocal

rem Editable parameters
set NAME="caenorhabditis elegans"
set WINDOW_LENGTH=300000
set DIR="caenorhabditis_elegans"
set GCF=GCF_000002985.6
set GTF_PATH=resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf
set PARTITIONS=300
set REGIONS=3
set K_RANGE="(4,12)"
set GENES_KMER_SIZE=4
set GENES_AMOUNT=100

rem Execute commands
echo Running download...
py .\command.py download -name %NAME%

echo Running analyze_genome -mode whole...
py .\command.py analyze_genome -name %NAME% -mode whole

echo Running analyze_genome -mode regions...
py .\command.py analyze_genome -name %NAME% -mode regions -window_length %WINDOW_LENGTH%

echo Running find_kmers_genome linearly...
py .\command.py find_kmers_genome -method l -k_range %K_RANGE% -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR% -graph_from_file true

echo Running graph -mode whole...
py .\command.py graph -name %NAME% -mode whole

echo Running graph -mode region...
py .\command.py graph -name %NAME% -mode region -window_length %WINDOW_LENGTH%

echo Running graph_linear_repeats_genome...
py .\command.py graph_linear_repeats_genome --save true -name %NAME% -gcf %GCF% -dir %DIR% -k_range %K_RANGE%

echo Running graph_gtf_file...
py .\command.py graph_gtf_file -path %GTF_PATH% -partitions %PARTITIONS% -regions %REGIONS% -plot_type line -dir %DIR% --save true

echo Running graph_linear_regression_genome...
py .\command.py graph_linear_regression_genome -gcf %GCF% -k_range %K_RANGE% -dir %DIR% -name %NAME%

echo Running load_genes...
py .\command.py load_genes -path %GTF_PATH%

echo Running find_kmers_linearly_genes_genome...
py .\command.py find_kmers_linearly_genes_genome -name %NAME% -window_length %WINDOW_LENGTH% -dir %DIR% -graph_from_file true -size %GENES_KMER_SIZE%

echo Scraping categories and subcategories..
py .\command.py load_categories_genome -name "caenorhabditis elegans" -size %GENES_KMER_SIZE% -genes_amount %GENES_AMOUNT%

:: Graph heatmaps

endlocal