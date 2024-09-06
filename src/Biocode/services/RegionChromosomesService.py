from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class RegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "region_chromosomes"
        self.columns = ["name", "refseq_accession_number", "organism_id", "cover_percentage", "cover", "regions_total", "region_number",
                        "size", "whole_chromosome_id"]
        self.pk_column = "id"

    def extract_id_by_refseq_accession_number(self, refseq_accession_number: str) -> int | None:
        result = self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)

        if result is not None and not result.empty:
            return int(result.loc[0, 'id'])
        return None