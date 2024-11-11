NAME="caenorhabditis remanei"
WINDOW_LENGTH=300000
DIR="caenorhabditis_remanei"
GCF=GCF_010183535.1
GTF_PATH=resources/genes/caenorhabditis_remanei/gtf/GCF_010183535.1_CRPX506_genomic.gtf
SEQUENCE_PATH=resources/dna_sequences/caenorhabditis_remanei/chrX.fna
PARTITIONS=300
REGIONS=3
K_RANGE="(4,12)"
GENES_KMER_SIZE=4
GENES_AMOUNT=100

.PHONY: all download analyze_sequence_whole analyze_sequence_regions \
        find_kmers_sequence graph_whole graph_region graph_linear_repeats \
        graph_gtf_file graph_linear_regression load_genes find_kmers_genes \
        graph_genes_repeats load_categories graph_categories

all: analyze_sequence_whole analyze_sequence_regions find_kmers_sequence \
     graph_whole graph_region graph_linear_repeats graph_gtf_file \
     graph_linear_regression load_genes find_kmers_genes graph_genes_repeats \
     load_categories graph_categories

download:
	@echo "Running download..."
	python3 command.py download -name $(NAME)

analyze_sequence_whole:
	@echo "Running analyze_sequence -mode whole..."
	python3 command.py analyze_sequence -path $(SEQUENCE_PATH) -name $(NAME) -mode whole

analyze_sequence_regions:
	@echo "Running analyze_sequence -mode regions..."
	python3 command.py analyze_sequence -path $(SEQUENCE_PATH) -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

find_kmers_sequence:
	@echo "Running find_kmers_sequence linearly..."
	python3 command.py find_kmers_sequence -path $(SEQUENCE_PATH) -method l -k_range $(K_RANGE) -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true

graph_whole:
	@echo "Running graph -mode whole..."
	python3 command.py graph -name $(NAME) -mode whole

graph_region:
	@echo "Running graph -mode region..."
	python3 command.py graph -name $(NAME) -mode region -window_length $(WINDOW_LENGTH)

graph_linear_repeats:
	@echo "Running graph_linear_repeats_sequence..."
	python3 command.py graph_linear_repeats_sequence --save true -name $(NAME) -path $(SEQUENCE_PATH) -dir $(DIR) -k_range $(K_RANGE)

graph_gtf_file:
	@echo "Running graph_gtf_file..."
	python3 command.py graph_gtf_file -path $(GTF_PATH) -partitions $(PARTITIONS) -regions $(REGIONS) -plot_type line -dir $(DIR) --save true

graph_linear_regression:
	@echo "Running graph_linear_regression_sequence..."
	python3 command.py graph_linear_regression_sequence -path $(SEQUENCE_PATH) -k_range $(K_RANGE) -dir $(DIR) -name $(NAME)

load_genes:
	@echo "Running load_genes..."
	python3 command.py load_genes -path $(GTF_PATH)

find_kmers_genes: load_genes
	@echo "Running find_kmers_linearly_genes_sequence..."
	python3 command.py find_kmers_linearly_genes_sequence -path $(SEQUENCE_PATH) -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true -size $(GENES_KMER_SIZE)

graph_genes_repeats: find_kmers_genes
	@echo "Running graph_linear_in_genes_repeats_sequence..."
	python3 command.py graph_linear_in_genes_repeats_sequence --save true -name $(NAME) -path $(SEQUENCE_PATH) -dir $(DIR) -k_range $(K_RANGE)

load_categories:
	@echo "Running load_categories..."
	python3 command.py load_categories -path $(SEQUENCE_PATH) -name $(NAME) -size $(GENES_KMER_SIZE) -genes_amount $(GENES_AMOUNT)

graph_categories: load_categories
	@echo "Running graph_categories_repeats_heatmap_sequence..."
	python3 command.py graph_categories_repeats_heatmap_sequence -path $(SEQUENCE_PATH) -size $(GENES_KMER_SIZE) --save true -dir $(DIR) -name $(NAME) --tags true
