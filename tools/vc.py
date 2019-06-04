import os
import loggers
logger = loggers.get('vc')
def git_init(path='.'):
    logger.info('Initializing git repo at %s...', path)
    os.system("git init %s" % path)
