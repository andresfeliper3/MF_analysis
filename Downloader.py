import os
import subprocess
import gzip
import shutil

from load import loader
from utils.decorators import Timer, TryExcept
from utils.logger import logger




class Downloader:

    def __init__(self):
        self.folder_path = "resources/dna_sequences/"
        self.organism = ""

    @TryExcept
    def download_command(self, args):
        if args.name:
            self.organism = args.name
            loader.set_organism(self.organism)
            self.remove_files(folder=loader.get_organism_folder())
            self.download_from_repository(folder=loader.get_organism_folder(),
                                                download_url=loader.get_download_url(),
                                                suffix=".gz")
            self.uncompress_all_files(folder=loader.get_organism_folder())
            if bool(args.gff):
                self.download_from_repository(folder=loader.get_organism_gtf_subfolder(),
                                                    download_url=loader.get_download_gff_url(), suffix=".gtf.gz")
                self.uncompress_all_files(folder=loader.get_organism_gtf_subfolder())
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

    @Timer
    def remove_files(self, folder):
        directory_path = os.path.join(self.folder_path, folder)

        if os.path.exists(directory_path) and os.path.isdir(directory_path):
            try:
                # Walk thr ough all the files and folders within the directory and delete them
                for root, dirs, files in os.walk(directory_path, topdown=False):
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
    def download_from_repository(self,folder: str, download_url: str, suffix: str):
        directory_path = f"{self.folder_path}{folder}"
        download_command = f'wget --recursive -np -e robots=off --reject "index.html" --no-host-directories ' \
                        f'--cut-dirs=10 --accept "*{suffix}" {download_url} -P {directory_path}'

        result = subprocess.run(download_command, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            logger.info(result.stdout)
        else:
            logger.error(result.stderr)


    @Timer
    def clean_directory(self, folder: str, suffix: str):
        directory_path = f"{self.folder_path}{folder}"
        files = os.listdir(directory_path)

        # Keep only files ending with suffix ending and delete others
        for file in files:
            file_path = os.path.join(directory_path, file)
            if not file.endswith(suffix):
                os.remove(file_path)


    @Timer
    def uncompress_all_files(self, folder):
        directory_path = f"{self.folder_path}{folder}"
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



