import os
import subprocess
import gzip
import shutil

from load import loader
from utils.decorators import Timer
from utils.logger import logger


FOLDER_PATH = "resources/dna_sequences/"

@Timer
def remove_files(folder):
    directory_path = os.path.join(FOLDER_PATH, folder)

    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        try:
            # Walk through all the files and folders within the directory
            for root, dirs, files in os.walk(directory_path, topdown=False):
                # Delete all files in the current directory
                for file in files:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    logger.info(f"Deleted file: {file_path}")

                # Delete all subdirectories in the current directory
                for dir in dirs:
                    dir_path = os.path.join(root, dir)
                    shutil.rmtree(dir_path)
                    logger.info(f"Deleted directory and its contents: {dir_path}")

            # Finally, delete the main directory itself if required
            shutil.rmtree(directory_path)
            logger.info(f"Deleted directory and its contents: {directory_path}")

        except Exception as e:
            logger.error(f"Error deleting contents of {directory_path}: {e}")
    else:
        logger.error(f"The directory {directory_path} does not exist.")


@Timer
def execute_download_command(folder: str, download_url: str, suffix: str):
    directory_path = f"{FOLDER_PATH}{folder}"


    download_command = f'wget --recursive -np -e robots=off --reject "index.html" --no-host-directories ' \
                    f'--cut-dirs=10 --accept "*{suffix}" {download_url} -P {directory_path}'

    result = subprocess.run(download_command, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        logger.info(result.stdout)
    else:
        logger.error(result.stderr)


@Timer
def clean_directory(folder: str, suffix: str):
    directory_path = f"{FOLDER_PATH}{folder}"
    files = os.listdir(directory_path)

    # Keep only files ending with suffix ending and delete others
    for file in files:
        file_path = os.path.join(directory_path, file)
        if not file.endswith(suffix):
            os.remove(file_path)


@Timer
def uncompress_all_files(folder):
    directory_path = f"{FOLDER_PATH}{folder}"
    files = os.listdir(directory_path)

    for file in files:
        if file.endswith(".gz"):
            input_path = os.path.join(directory_path, file)
            output_path = os.path.join(directory_path, file.replace(".gz", ""))

            with gzip.open(input_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())

            # Optionally, you can delete the compressed file if needed
            os.remove(input_path)
            logger.info(f"Uncompressed file: {output_path}")



