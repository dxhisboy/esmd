import logging
import uuid
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
        levelname = record.levelname
        record.levelname = COLOR_LEVEL_NAME[record.levelname]
        ret = super().format(record)
        record.levelname = levelname
        return ret

fformat = logging.Formatter('%(asctime)s - %(name)s:%(levelname)s - %(message)s')
cformat = ColoredFormatter('%(asctime)s - %(name)s:%(levelname)s - %(message)s')

def add_file(name, logfile):
    logger = logging.getLogger(None)
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fformat)
    logger.addHandler(fh)

def init_root_logger():
    logger = logging.getLogger(None)
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(cformat)
    logger.addHandler(ch)

def get(name = None):
    return logging.getLogger(name)

if __name__ == 'loggers':
    init_root_logger()
