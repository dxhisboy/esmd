import logging
import time
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(30, 38)

COLOR_SEQ = "\033[1;%dm%s\033[0m"

COLOR_LEVEL_NAME = {
    'WARNING': COLOR_SEQ % (YELLOW, 'WARNING'),
    'INFO': COLOR_SEQ % (GREEN, 'INFO'),
    'DEBUG': COLOR_SEQ % (CYAN, 'DEBUG'),
    'CRITICAL': COLOR_SEQ % (MAGENTA, 'CRITICAL'),
    'ERROR': COLOR_SEQ % (RED, 'ERROR'),
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%'):
        super().__init__(fmt, datefmt)
    def format(self, record):
        record.levelname = COLOR_LEVEL_NAME[record.levelname]
        return super().format(record)

fformat = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
cformat = ColoredFormatter('%(asctime)s - %(levelname)s - %(message)s')

def create(logfile = None):
    logger = logging.getLogger('simple_example')
    logger.setLevel(logging.INFO)
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fformat)
        logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(cformat)
    logger.addHandler(ch)

    return logger

logger = create('test.log')
# 'application' code
logger.debug('debug message')
logger.info('info message')
logger.warning('warn message')
logger.error('error message')
logger.critical('critical message')
