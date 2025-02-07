@echo off
setlocal

rem Editable parameters
set NAME="caenorhabditis elegans"
set WINDOW_LENGTH=300000
set DIR="caenorhabditis_elegans"
set GCF=GCF_000002985.6
set PARTITIONS=100
set REGIONS=3
set K_RANGE="(4,12)"
set GENES_KMER_SIZE=4
set GENES_AMOUNT=100

echo Running graph_cgr -mode whole...
py .\command.py graph_cgr -name %NAME% -mode whole

::echo Running analyze_genome -mode regions...
::py .\command.py analyze_genome -name %NAME% -mode regions -window_length %WINDOW_LENGTH%


endlocal