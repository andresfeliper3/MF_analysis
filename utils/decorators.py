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

def Singleton(cls):
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            logger.info(f"Creating instance of {cls.__name__}")
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]

    return get_instance

def Inject(**dependencies):
    def decorator(cls):
        original_init = cls.__init__

        def __init__(self, *args, **kwargs):
            for key, dependency in dependencies.items():
                if key not in kwargs:
                    kwargs[key] = dependency()
            original_init(self, *args, **kwargs)

        cls.__init__ = __init__
        return cls

    return decorator