NAME="aspergillus luchuensis"
WINDOW_LENGTH=150000
DIR="aspergillus_luchuensis"
GCF=GCF_016861625.1
GTF_PATH=resources/genes/aspergillus_luchuensis/gtf/GCF_016861625.1_AkawachiiIFO4308_assembly01_genomic.gtf
PARTITIONS=100
REGIONS=3
K_RANGE="(4,12)"
GENES_KMER_SIZE=4

# Define targets
.PHONY: all graph_whole graph_regions graph_linear_repeats graph_gtf_file graph_linear_regression graph_categories

all: graph_whole graph_regions graph_linear_repeats graph_gtf_file graph_linear_regression graph_categories

graph_whole:
	@echo "Running graph -mode whole..."
	python3 command.py graph -name $(NAME) -mode whole

graph_regions:
	@echo "Running graph -mode regions..."
	python3 command.py graph -name $(NAME) -mode regions -window_length $(WINDOW_LENGTH)

graph_linear_repeats:
	@echo "Running graph_linear_repeats_genome..."
	python3 command.py graph_linear_repeats_genome --save true -name $(NAME) -gcf $(GCF) -dir $(DIR) -k_range $(K_RANGE)

graph_gtf_file:
	@echo "Running graph_gtf_file..."
	python3 command.py graph_gtf_file -path $(GTF_PATH) -partitions $(PARTITIONS) -regions $(REGIONS) -plot_type line -dir $(DIR) --save true

graph_linear_regression:
	@echo "Running graph_linear_regression_genome..."
	python3 command.py graph_linear_regression_genome -gcf $(GCF) -k_range $(K_RANGE) -dir $(DIR) -name $(NAME)

graph_categories:
	@echo "Running graph_categories_repeats_heatmap_genome..."
	python3 command.py graph_categories_repeats_heatmap_genome -size $(GENES_KMER_SIZE) --save true -dir $(DIR) -name $(NAME) --tags true
