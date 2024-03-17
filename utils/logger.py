import logging

# Configure the logging module
logging.basicConfig(level=logging.DEBUG,  # Set the minimum log level
                    format='%(asctime)s - %(levelname)s - %(message)s')  # Define the log message format

# Create a logger
logger = logging.getLogger(__name__)

