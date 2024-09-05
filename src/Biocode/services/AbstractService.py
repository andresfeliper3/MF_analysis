from src.Biocode.managers.DBConnectionManager import DBConnectionManager


class AbstractService:

    def insert(self, record: tuple) -> int :
        return DBConnectionManager.insert(table_name=self.table_name, columns=self.columns, pk=self.pk_column,
                                          record=record)
    def update(self, pk_value, record: tuple) -> int:
        return DBConnectionManager.update_when_null(table_name=self.table_name, columns=self.columns, pk_col=self.pk_column,
                                          pk_value=pk_value, record=record)
    def extract_all(self):
        return DBConnectionManager.extract_all(table_name=self.table_name)

    def extract_by_id(self, target_id: int):
        return DBConnectionManager.extract_by_target(table_name=self.table_name, column=self.pk_column,
                                                     target=target_id)

    def extract_by_field(self, column: str, value):
        return DBConnectionManager.extract_by_target(table_name=self.table_name, column=column, target=value)

    def extract_with_custom_query(self, query):
        return DBConnectionManager.extract_with_custom_query(query)