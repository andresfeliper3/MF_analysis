NAME="malassezia restricta"
WINDOW_LENGTH=30000
DIR="malassezia_restricta"
GCF=GCF_003290485.1
GTF_PATH=resources/genes/malassezia_restricta/gtf/GCF_003290485.1_ASM329048v1_genomic.gtf
PARTITIONS=300
REGIONS=3
K_RANGE="(4,12)"
GENES_KMER_SIZE=4
GENES_AMOUNT=100

.PHONY: all download analyze_genome_whole analyze_genome_regions \
        graph_whole graph_regions find_kmers_genome \
        graph_linear_repeats graph_gtf_file graph_linear_regression \
        load_genes find_kmers_genes load_categories graph_categories

all: analyze_genome_whole analyze_genome_regions graph_whole \
     graph_regions find_kmers_genome graph_linear_repeats \
     graph_gtf_file graph_linear_regression load_genes \
     find_kmers_genes load_categories graph_categories

download:
	@echo "Running download..."
	python3 command.py download -name $(NAME)

analyze_genome_whole:
	@echo "Running analyze_genome -mode whole..."
	python3 command.py analyze_genome -name $(NAME) -mode whole

analyze_genome_regions:
	@echo "Running analyze_genome -mode regions..."
	python3 command.py analyze_genome -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

graph_whole:
	@echo "Running graph -mode whole..."
	python3 command.py graph -name $(NAME) -mode whole

graph_regions:
	@echo "Running graph -mode regions..."
	python3 command.py graph -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

find_kmers_genome:
	@echo "Running find_kmers_genome linearly..."
	python3 command.py find_kmers_genome -method l -k_range $(K_RANGE) -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR)

graph_linear_repeats:
	@echo "Running graph_linear_repeats_genome..."
	python3 command.py graph_linear_repeats_genome --save true -name $(NAME) -gcf $(GCF) -dir $(DIR) -k_range $(K_RANGE)

graph_gtf_file:
	@echo "Running graph_gtf_file..."
	python3 command.py graph_gtf_file -path $(GTF_PATH) -partitions $(PARTITIONS) -regions $(REGIONS) -plot_type line -dir $(DIR) --save true

graph_linear_regression:
	@echo "Running graph_linear_regression_genome..."
	python3 command.py graph_linear_regression_genome -gcf $(GCF) -k_range $(K_RANGE) -dir $(DIR) -name $(NAME)

load_genes:
	@echo "Running load_genes..."
	python3 command.py load_genes -path $(GTF_PATH)

find_kmers_genes: load_genes
	@echo "Running find_kmers_linearly_genes_genome..."
	python3 command.py find_kmers_linearly_genes_genome -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true -size $(GENES_KMER_SIZE)

load_categories:
	@echo "Running load_categories_genome..."
	python3 command.py load_categories_genome -name $(NAME) -size $(GENES_KMER_SIZE) -genes_amount $(GENES_AMOUNT)

graph_categories: load_categories
	@echo "Running graph_categories_repeats_heatmap_genome..."
	python3 command.py graph_categories_repeats_heatmap_genome -size $(GENES_KMER_SIZE) --save true -dir $(DIR) -name $(NAME) --tags true
