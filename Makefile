.PHONY: all clean

# Editable parameters
NAME ?= "caenorhabditis elegans"
WINDOW_LENGTH ?= 300000
DIR ?= "caenorhabditis_elegans"
GCF ?= GCF_000002985.6
GTF_PATH ?= resources/genes/caenorhabditis_elegans/gtf/GCF_000002985.6_WBcel235_genomic.gtf
PARTITIONS ?= 300
REGIONS ?= 3
K_RANGE ?= "(4,12)"
GENES_KMER_SIZE ?= 4
GENES_AMOUNT ?= 100

all: download analyze_genome_whole analyze_genome_regions find_kmers_genome graph_whole graph_region graph_linear_repeats_genome graph_gtf_file graph_linear_regression_genome load_genes find_kmers_linearly_genes_genome load_categories_genome graph_categories_repeats_heatmap_genome

download:
	@echo "Running download..."
	py ./command.py download -name $(NAME)

analyze_genome_whole:
	@echo "Running analyze_genome -mode whole..."
	py ./command.py analyze_genome -name $(NAME) -mode whole

analyze_genome_regions:
	@echo "Running analyze_genome -mode regions..."
	py ./command.py analyze_genome -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

find_kmers_genome:
	@echo "Running find_kmers_genome..."
	py ./command.py find_kmers_genome -method l -k_range $(K_RANGE) -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true

graph_whole:
	@echo "Running graph -mode whole..."
	py ./command.py graph -name $(NAME) -mode whole

graph_region:
	@echo "Running graph -mode region..."
	py ./command.py graph -name $(NAME) -mode region -window_length $(WINDOW_LENGTH)

graph_linear_repeats_genome:
	@echo "Running graph_linear_repeats_genome..."
	py ./command.py graph_linear_repeats_genome --save true -name $(NAME) -gcf $(GCF) -dir $(DIR) -k_range $(K_RANGE)

graph_gtf_file:
	@echo "Running graph_gtf_file..."
	py ./command.py graph_gtf_file -path $(GTF_PATH) -partitions $(PARTITIONS) -regions $(REGIONS) -plot_type line -dir $(DIR) --save true

graph_linear_regression_genome:
	@echo "Running graph_linear_regression_genome..."
	py ./command.py graph_linear_regression_genome -gcf $(GCF) -k_range $(K_RANGE) -dir $(DIR) -name $(NAME)

load_genes:
	@echo "Running load_genes..."
	py ./command.py load_genes -path $(GTF_PATH)

find_kmers_linearly_genes_genome:
	@echo "Running find_kmers_linearly_genes_genome..."
	py ./command.py find_kmers_linearly_genes_genome -name $(NAME) -window_length $(WINDOW_LENGTH) -dir $(DIR) -graph_from_file true -size $(GENES_KMER_SIZE)

load_categories_genome:
	@echo "Scraping functional categories and subcategories from KEGG..."
	py ./command.py load_categories_genome -name $(NAME) -size $(GENES_KMER_SIZE) -genes_amount $(GENES_AMOUNT)

graph_categories_repeats_heatmap_genome:
	@echo "Graphing functional categories and subcategories heatmaps..."
	py ./command.py graph_categories_repeats_heatmap_genome -size $(GENES_KMER_SIZE) --save true -dir $(DIR) -name $(NAME) --tags true

clean:
	@echo "Cleaning up..."
	# Add commands to clean up generated files