from utils.decorators import Singleton
from utils.logger import logger


def Service(cls):
    singleton_cls = Singleton(cls)

    # Add additional features if necessary
    def get_instance(*args, **kwargs):
        return singleton_cls(*args, **kwargs)

    return get_instance





