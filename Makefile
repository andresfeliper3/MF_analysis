.PHONY: all clean

# Makefile for running genome analysis tasks

# Editable parameters
NAME := "drosophila pseudoobscura"
WINDOW_LENGTH := 600000
DIR := "drosophila_pseudoobscura"
GCF := GCF_009870125.1
GTF_PATH := resources/genes/drosophila_pseudoobscura/gtf/GCF_009870125.1_UCI_Dpse_MV25_genomic.gtf
PARTITIONS := 300
REGIONS := 3
K_RANGE := "(4,12)"
GENES_KMER_SIZE := 4
GENES_AMOUNT := 100

.PHONY: all analyze_genome find_kmers_genome graph load_genes find_kmers_linearly_genes_genome load_categories_genome graph_categories_repeats_heatmap_genome

all: analyze_genome find_kmers_genome graph load_genes find_kmers_linearly_genes_genome load_categories_genome graph_categories_repeats_heatmap_genome

analyze_genome:
	@echo "Running analyze_genome -mode whole..."
	py ./command.py analyze_genome -name $(NAME) -mode whole
	@echo "Running analyze_genome -mode regions..."
	py ./command.py analyze_genome -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

find_kmers_genome:
	@echo "Running find_kmers_genome linearly..."
	py ./command.py find_kmers_genome -method l -k_range $(K_RANGE) -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR)

graph:
	@echo "Running graph -mode whole..."
	py ./command.py graph -name $(NAME) -mode whole
	@echo "Running graph -mode regions..."
	py ./command.py graph -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)
	@echo "Running graph_linear_repeats_genome..."
	py ./command.py graph_linear_repeats_genome --save true -name $(NAME) -gcf $(GCF) -dir $(DIR) -k_range $(K_RANGE)
	@echo "Running graph_gtf_file..."
	py ./command.py graph_gtf_file -path $(GTF_PATH) -partitions $(PARTITIONS) -regions $(REGIONS) -plot_type line -dir $(DIR) --save true
	@echo "Running graph_linear_regression_genome..."
	py ./command.py graph_linear_regression_genome -gcf $(GCF) -k_range $(K_RANGE) -dir $(DIR) -name $(NAME)

load_genes:
	@echo "Running load_genes..."
	py ./command.py load_genes -path $(GTF_PATH)

find_kmers_linearly_genes_genome:
	@echo "Running find_kmers_linearly_genes_genome..."
	py ./command.py find_kmers_linearly_genes_genome -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true -size $(GENES_KMER_SIZE)

load_categories_genome:
	@echo "Running scraping functional categories and subcategories from KEGG..."
	py ./command.py load_categories_genome -name $(NAME) -size $(GENES_KMER_SIZE) -genes_amount $(GENES_AMOUNT)

graph_categories_repeats_heatmap_genome:
	@echo "Running Graphing functional categories and subcategories heatmaps..."
	py ./command.py graph_categories_repeats_heatmap_genome -size $(GENES_KMER_SIZE) --save true -dir $(DIR) -name $(NAME) --tags true

clean:
	@echo "Cleaning up..."
	# Add commands to clean up generated files