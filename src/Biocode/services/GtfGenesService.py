from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GtfGenesService(AbstractService):
    def __init__(self):
        self.table_name = "gtf_genes"
        self.columns = ["whole_chromosomes_id", "source", "feature", "start_position", "end_position", "length",
                        "score", "strand", "frame", "gene_id_gtf", "gene", "gene_biotype"]
        self.pk_column = "id"

    def extract_genes_by_chromosome(self, refseq_accession_number: str):
        query = f"SELECT gtf_genes.id, source, feature, start_position, end_position, length, score, strand, frame, " \
                f"gene_id_gtf, gene, gene_biotype, refseq_accession_number, size, o.name FROM gtf_genes JOIN " \
                f"whole_chromosomes wc on gtf_genes.whole_chromosomes_id = wc.id JOIN organisms o " \
                f"on wc.organism_id = o.id WHERE refseq_accession_number='{refseq_accession_number}';"
        return self.extract_with_custom_query(query)

