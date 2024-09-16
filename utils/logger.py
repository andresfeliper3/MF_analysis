# logging_config.py
import logging
import os
from datetime import datetime

# Define the folder for log files
log_folder = "logs"

# Ensure the log folder exists
if not os.path.exists(log_folder):
    os.makedirs(log_folder)

# Generate a log file name with the current date
log_filename = datetime.now().strftime("%Y-%m-%d") + ".log"
log_file_path = os.path.join(log_folder, log_filename)

# Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set the minimum log level

# Prevent duplicate handlers if this module is imported multiple times
if not logger.hasHandlers():
    # Create handlers
    console_handler = logging.StreamHandler()  # Handler for console output
    file_handler = logging.FileHandler(log_file_path)  # Handler for log file with date-based file name

    # Set levels for handlers
    console_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.INFO)

    # Create formatters
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Assign formatters to handlers
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
