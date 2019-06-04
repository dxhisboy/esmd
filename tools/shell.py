import subprocess
import loggers
logger = loggers.get('shell')
def run(cmd, shell=True):
    logger.debug("executing: %s" % cmd)
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=shell)
    logger.debug("========command output========\n%s" % output.decode().rstrip())
    logger.debug("==========end output==========")
    
