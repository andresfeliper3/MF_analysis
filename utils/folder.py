import os
from utils.logger import logger


def apply_function_to_files_in_folder(directory_path: str, func, *args, **kwargs):
    """
    Apply a given function to each file in a specified directory.

    :param directory_path: Path to the directory containing the files.
    :param func: The function to apply to each file.
    :param args: Positional arguments to pass to the function.
    :param kwargs: Keyword arguments to pass to the function.
    """
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        files = os.listdir(directory_path)
        for file in files:
            file_path = os.path.join(directory_path, file)
            func(file_path, *args, **kwargs)
    else:
        logger.error(f"The directory {directory_path} does not exist.")
