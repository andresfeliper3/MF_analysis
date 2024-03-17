import pandas as pd
import os
import sqlite3


class DBConnectionManager:
    db_directory = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')),
                                "MFA_analysis.db")

    @staticmethod
    def start():
        DBConnectionManager.conn = sqlite3.connect(DBConnectionManager.db_directory)
        DBConnectionManager.cursor = DBConnectionManager.conn.cursor()

    """
    Example usage:
    table_name = "organisms"
    columns = ["name", "GCF"]
    data = [("Organism1", "123"), ("Organism2", "456")]
    """

    @staticmethod
    def insert(table_name, columns: list, pk, record: tuple):
        check_query = f'SELECT * FROM {table_name} WHERE {" AND ".join([f"{col} = ?" for col in columns if col != pk])};'
        DBConnectionManager.cursor.execute(check_query, tuple(record))
        existing_record = DBConnectionManager.cursor.fetchone()

        if existing_record:
            existing_record_id = existing_record[0]  # Assuming the id is the first element of the tuple
            print(f"Record {record} already exists. Skipping insertion.")
            print("existing record")
            print(existing_record)
            return existing_record_id
        else:
            insert_query = f'INSERT INTO {table_name} ({", ".join(columns)}) VALUES ({", ".join(["?"] * len(columns))});'
            DBConnectionManager.cursor.execute(insert_query, tuple(record))

        DBConnectionManager.conn.commit()
        last_inserted_id = DBConnectionManager.cursor.lastrowid
        return last_inserted_id

    @staticmethod
    def extract_all(table_name):
        query = f"SELECT * FROM {table_name};"
        DBConnectionManager.cursor.execute(query)
        rows = DBConnectionManager.cursor.fetchall()

        if not rows:
            return None

        columns = [col[0] for col in DBConnectionManager.cursor.description]
        df = pd.DataFrame(rows, columns=columns)
        return df

    @staticmethod
    def extract_by_target(table_name, column, target):
        query = f"SELECT * FROM {table_name} WHERE {column} = ?;"
        DBConnectionManager.cursor.execute(query, (target,))
        rows =DBConnectionManager.cursor.fetchall()

        if not rows:
            return None

        columns = [col[0] for col in DBConnectionManager.cursor.description]
        df = pd.DataFrame(rows, columns=columns)
        return df

    @staticmethod
    def extract_with_custom_query(query):
        DBConnectionManager.cursor.execute(query)
        rows = DBConnectionManager.cursor.fetchall()

        if not rows:
            return None

        columns = [col[0] for col in DBConnectionManager.cursor.description]
        df = pd.DataFrame(rows, columns=columns)
        return df

    @staticmethod
    def close():
        DBConnectionManager.conn.close()
