from src.Biocode.services.AbstractService import AbstractService

class GtfGenesService(AbstractService):
    def __init__(self):
        self.table_name = "gtf_genes"
        self.columns = ["whole_chromosomes_id", "source", "feature", "start_position", "end_position", "length",
                        "score", "strand", "frame", "gene_id_gtf", "gene", "gene_biotype"]
        self.pk_column = "id"

    def extract_genes_by_genome(self, GCF: str):
        self.ex