from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class WholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "whole_chromosomes"
        self.columns = ["name", "refseq_accession_number", "organism_id", "cover_percentage", "cover", "size"]
        self.pk_column = "id"


    def extract_by_name(self, name: str) -> str:
        return self.extract_by_field(column="name", value=name)

    def extract_by_refseq_accession_number(self, refseq_accession_number: str) -> str:
        return self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)

    def extract_id_by_refseq_accession_number(self, refseq_accession_number: str) -> int:
        return int(self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number).loc[0, 'id'])

