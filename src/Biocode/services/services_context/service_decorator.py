from utils.decorators import Singleton
from utils.logger import logger


def Service(cls):
    singleton_cls = Singleton(cls)

    # Add additional features if necessary
    cls.get_instance = lambda *args, **kwargs: singleton_cls(*args, **kwargs)

    return cls






