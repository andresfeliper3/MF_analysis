from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

from typing import List

@Service
class OrganismsService(AbstractService):

    def __init__(self):
        self.table_name = "organisms"
        self.columns = ["name", "GCF", "amount_chromosomes"]
        self.pk_column = "id"

    def extract_by_GCF(self, GCF: str):
        return int(self.extract_by_field(column="GCF", value=GCF).loc[0, 'id'])

    def extract_chromosomes_refseq_accession_numbers_by_GCF(self, GCF: str) -> List:
        query = f"SELECT refseq_accession_number FROM whole_chromosomes JOIN organisms o on whole_chromosomes.organism_id = o.id " \
                f"WHERE GCF='{GCF}';"
        return self.extract_with_custom_query(query)['refseq_accession_number'].to_numpy()

    def extract_chromosomes_names_by_GCF(self, GCF: str) -> List:
        query = f"SELECT whole_chromosomes.name FROM whole_chromosomes JOIN organisms o on whole_chromosomes.organism_id = o.id " \
                f"WHERE GCF='{GCF}';"
        return self.extract_with_custom_query(query)['name'].to_numpy()

