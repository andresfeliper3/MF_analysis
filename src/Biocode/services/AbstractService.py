from src.Biocode.managers.DBConnectionManager import DBConnectionManager


class AbstractService:

    def insert(self, record: tuple):
        return DBConnectionManager.insert(table_name=self.table_name, columns=self.columns, pk=self.pk_column,
                                          record=record)

    def extract_all(self):
        return DBConnectionManager.extract_all(table_name=self.table_name)

    def extract_by_id(self, target_id: int):
        return DBConnectionManager.extract_by_target(table_name=self.table_name, column=self.pk_column,
                                                     target=target_id)
