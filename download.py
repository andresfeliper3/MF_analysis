import os
import subprocess
import gzip

from load import loader
from utils.decorators import Timer
from utils.logger import logger


FOLDER_PATH = "resources/dna_sequences/"

@Timer
def remove_files(organism_folder):
    directory_path = f"{FOLDER_PATH}{organism_folder}"

    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        files = os.listdir(directory_path)

        for file in files:
            file_path = os.path.join(directory_path, file)

            try:
                os.remove(file_path)
                logger.info(f"Deleted: {file_path}")
            except Exception as e:
                logger.error(f"Error deleting {file_path}: {e}")
    else:
        logger.error(f"The directory {directory_path} does not exist.")


@Timer
def execute_download_command(organism_folder, download_url):
    directory_path = f"{FOLDER_PATH}{organism_folder}"


    download_command = f'wget --recursive -np -e robots=off --reject "index.html" --no-host-directories ' \
                    f'--cut-dirs=10 {download_url} -P {directory_path}'

    result = subprocess.run(download_command, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        logger.info(result.stdout)
    else:
        logger.error(result.stderr)


@Timer
def clean_directory(organism_folder):
    directory_path = f"{FOLDER_PATH}{organism_folder}"
    files = os.listdir(directory_path)

    # Keep only files ending with ".gz" and delete others
    for file in files:
        file_path = os.path.join(directory_path, file)
        if not file.endswith(".gz"):
            os.remove(file_path)


@Timer
def uncompress_all_files(organism_folder):
    directory_path = f"{FOLDER_PATH}{organism_folder}"
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



