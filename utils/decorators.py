import os
import sys
from src.Biocode.managers.DBConnectionManager import DBConnectionManager

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import time
from utils.logger import logger


def Timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"Execution time for {func.__name__}: {execution_time} seconds")
        return result
    return wrapper

def DBConnection(func):
    def wrapper(*args, **kwargs):
        DBConnectionManager.start()
        logger.info(f"Database connection started")
        func(*args, **kwargs)
        DBConnectionManager.close()
        logger.info(f"Database connection closed")
    return wrapper