import logging

# Configure the logging module
logging.basicConfig(level=logging.INFO,  # Set the minimum log level
                    format='%(asctime)s - %(levelname)s - %(message)s')  # Define the log message format

logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Create a logger
logger = logging.getLogger(__name__)

