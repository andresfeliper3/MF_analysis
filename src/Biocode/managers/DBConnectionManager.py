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
    def _get_placeholder_and_conversion():
        """Get the SQL placeholder and field conversion logic based on the database type."""
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

        return placeholder, array_conversion

    @staticmethod
    def _execute_query(query: str, record: tuple, return_last_id: bool = False, pk: str = None):
        """Executes the provided query with the given record and returns the last inserted ID if requested."""
        try:
            DBConnectionManager.cursor.execute(query, record)
            DBConnectionManager.conn.commit()

            if return_last_id:
                if DBConnectionManager.is_postgresql:
                    return DBConnectionManager.cursor.fetchone()[0]
                elif DBConnectionManager.is_sqlite or DBConnectionManager.is_mysql:
                    return DBConnectionManager.cursor.lastrowid

        except Exception as e:
            logger.error(f"Error executing query: {e}")
            raise

    @staticmethod
    def insert(table_name, columns: list, pk, record: tuple):
        # Get placeholder and conversion logic
        placeholder, array_conversion = DBConnectionManager._get_placeholder_and_conversion()

        # Convert the record fields (handling lists/arrays if necessary)
        converted_record = tuple(array_conversion(field) for field in record)

        # Build and execute check query to find existing record
        check_query = f'SELECT * FROM {table_name} WHERE {" AND ".join([f"{col} = {placeholder}" for col in columns if col != pk])};'
        DBConnectionManager.cursor.execute(check_query, converted_record)
        existing_record = DBConnectionManager.cursor.fetchone()

        if existing_record:
            existing_record_id = existing_record[0]
            logger.info(f"Existing record with id {existing_record_id} in table {table_name}")
            return existing_record_id

        # Build the insert query
        insert_query = f'INSERT INTO {table_name} ({", ".join(columns)}) VALUES ({", ".join([placeholder] * len(columns))})'
        if DBConnectionManager.is_postgresql:
            insert_query += f' RETURNING {pk};'

        # Execute the insert query and return the last inserted ID
        return DBConnectionManager._execute_query(insert_query, converted_record, return_last_id=True, pk=pk)

    @staticmethod
    def update(table_name, columns: list, pk_col, pk_value, record: tuple):
        placeholder, array_conversion = DBConnectionManager._get_placeholder_and_conversion()
        converted_record = tuple(array_conversion(field) for field in record)
        update_query = f'UPDATE {table_name} SET {", ".join([f"{col} = {placeholder}" for col in columns if col != pk_col])} WHERE {pk_col} = {pk_value};'
        logger.info(f"Updating row in table {table_name} with {pk_col} = {pk_value}")
        return DBConnectionManager._execute_query(update_query, converted_record)


    @staticmethod
    def update_when_null(table_name, columns: list, pk_col, pk_value, record: tuple):
        placeholder, array_conversion = DBConnectionManager._get_placeholder_and_conversion()
        converted_record = tuple(array_conversion(field) for field in record)
        null_check_query = f"SELECT * FROM {table_name} WHERE {pk_col} = {placeholder} AND ({' OR '.join([f'{col} IS NULL' for col in columns])});"
        DBConnectionManager._execute_query(null_check_query, (pk_value,))
        row_with_nulls = DBConnectionManager.cursor.fetchone()  # Fetch one row with NULL values

        # If there is at least one row with a NULL value, proceed with the update
        if row_with_nulls:
            update_query = f'UPDATE {table_name} SET {", ".join([f"{col} = {placeholder}" for col in columns if col != pk_col])} WHERE {pk_col} = {pk_value};'
            logger.info(f"Updating row in table {table_name} with {pk_col} = {pk_value} where NULL values existed")
            return DBConnectionManager._execute_query(update_query, converted_record)
        else:
            logger.info(
                f"No NULL values found for row with {pk_col} = {pk_value} in table {table_name}. Update not performed.")

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
        try:
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
        except Exception as e:
            logger.error(f"Error executing query for {target} in {table_name}: {e}")
            raise

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
