import os
import yaml
import pandas as pd
from utils.logger import logger

class DBConnectionManager:
    conn = None
    cursor = None
    yaml_file = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')),
                                "resources/db_config.yaml")
    is_sqlite = False
    is_postgresql = False
    is_mysql = False

    @staticmethod
    def load_config(yaml_file):
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)

    @staticmethod
    def start():
        config = DBConnectionManager.load_config(DBConnectionManager.yaml_file)
        db_type = config['database']['type']

        if db_type == 'sqlite':
            import sqlite3
            db_file = config['database']['sqlite']['db_file']
            DBConnectionManager.conn = sqlite3.connect(db_file)
            DBConnectionManager.is_sqlite = True
        elif db_type == 'postgresql':
            import psycopg2
            DBConnectionManager.conn = psycopg2.connect(
                host=config['database']['postgresql']['host'],
                port=config['database']['postgresql']['port'],
                user=config['database']['postgresql']['user'],
                password=config['database']['postgresql']['password'],
                database=config['database']['postgresql']['database']
            )
            DBConnectionManager.is_postgresql = True
        elif db_type == 'mysql':
            import mysql.connector
            DBConnectionManager.conn = mysql.connector.connect(
                host=config['database']['mysql']['host'],
                port=config['database']['mysql']['port'],
                user=config['database']['mysql']['user'],
                password=config['database']['mysql']['password'],
                database=config['database']['mysql']['database']
            )
            DBConnectionManager.is_mysql = True
        else:
            raise ValueError(f"Unsupported database type: {db_type}")

        DBConnectionManager.cursor = DBConnectionManager.conn.cursor()

    @staticmethod
    def is_connected():
        return DBConnectionManager.conn is not None

    @staticmethod
    def insert(table_name, columns: list, pk, record: tuple):
        if DBConnectionManager.is_sqlite:
            placeholder = "?"
            array_conversion = lambda x: x if not isinstance(x, list) else str(x)
        elif DBConnectionManager.is_postgresql:
            placeholder = "%s"
            array_conversion = lambda x: x if not isinstance(x, list) else '{' + ','.join(map(str, x)) + '}'
        elif DBConnectionManager.is_mysql:
            placeholder = "%s"
            array_conversion = lambda x: x if not isinstance(x, list) else str(x)
        else:
            raise ValueError("Unsupported database type")

        converted_record = tuple(array_conversion(field) for field in record)

        check_query = f'SELECT * FROM {table_name} WHERE {" AND ".join([f"{col} = {placeholder}" for col in columns if col != pk])};'


        try:
            DBConnectionManager.cursor.execute(check_query, converted_record)
            existing_record = DBConnectionManager.cursor.fetchone()

            if existing_record:
                existing_record_id = existing_record[0]
                logger.info(f"Existing record with id {existing_record_id} in table {table_name}")
                return existing_record_id
            else:
                insert_query = f'INSERT INTO {table_name} ({", ".join(columns)}) VALUES ({", ".join([placeholder] * len(columns))})'

                if DBConnectionManager.is_postgresql:
                    insert_query += f' RETURNING {pk};'

                DBConnectionManager.cursor.execute(insert_query, converted_record)
                DBConnectionManager.conn.commit()
                logger.info(f"Inserted in table {table_name} - row {converted_record}")
                if DBConnectionManager.is_postgresql:
                    last_inserted_id = DBConnectionManager.cursor.fetchone()[0]
                elif DBConnectionManager.is_sqlite or DBConnectionManager.is_mysql:
                    last_inserted_id = DBConnectionManager.cursor.lastrowid

                return last_inserted_id
        except Exception as e:
            logger.error(f"Error executing query: {e}")
            raise

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
        # Determine the correct placeholder for the current database
        if DBConnectionManager.is_sqlite:
            placeholder = "?"
        else:
            placeholder = "%s"

        query = f"SELECT * FROM {table_name} WHERE {column} = {placeholder};"
        DBConnectionManager.cursor.execute(query, (target,))
        rows = DBConnectionManager.cursor.fetchall()

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
        if DBConnectionManager.conn is not None:
            DBConnectionManager.conn.close()
            DBConnectionManager.conn = None
            DBConnectionManager.cursor = None
        else:
            logger.warning("Database connection is not active.")
