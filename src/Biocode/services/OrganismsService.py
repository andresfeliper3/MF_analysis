from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.AbstractService import AbstractService


class OrganismsService(AbstractService):

    def __init__(self):
        self.table_name = "organisms"
        self.columns = ["name", "GCF", "amount_chromosomes"]
        self.pk_column = "id"

    def extract_by_GCF(self, GCF: str):
        return DBConnectionManager.extract_by_target(table_name=self.table_name, column="GCF", target=GCF)
