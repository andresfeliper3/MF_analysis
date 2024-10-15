from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

from typing import Tuple


@Service
class WholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "whole_chromosomes"
        self.columns = ["name", "refseq_accession_number", "organism_id", "cover_percentage", "cover", "size"]
        self.pk_column = "id"

    def extract_by_name(self, name: str) -> str | None:
        return self.extract_by_field(column="name", value=name)

    def extract_sequence_name_by_refseq_accession_number(self, refseq_accession_number: str):
        result = self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)

        if result is not None and not result.empty:
            return result.loc[0, 'name']
        return None

    def extract_by_refseq_accession_number(self, refseq_accession_number: str) -> str:
        return self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)

    def extract_id_by_refseq_accession_number(self, refseq_accession_number: str) -> int | None:
        result = self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)

        if result is not None and not result.empty:
            return int(result.loc[0, 'id'])
        return None

    def extract_size_by_refseq_accession_number(self, refseq_accession_number: str) -> int:
        return int(
            self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number).loc[0, 'size'])

    def extract_filename_by_refseq_accession_number(self, refseq_accession_number: str) -> str:
        return self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number).loc[0, 'name']

    def extract_filename_and_size_by_refseq_accession_number(self, refseq_accession_number: str) -> Tuple[str, int]:
        row = self.extract_by_field(column="refseq_accession_number", value=refseq_accession_number)
        return row.loc[0, 'name'], int(row.loc[0, 'size'])
