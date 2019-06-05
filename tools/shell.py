import subprocess
import logging
import loggers
import os
import threading

logger = loggers.get('shell')

class LoggerPipe(threading.Thread):
    def __init__(self, logger, level = logging.INFO):
        super().__init__()
        self.logger = logger
        self.level = level
        self.pipe = os.pipe()
        self.pin = os.fdopen(self.pipe[0])
        self.start()
    def fileno(self):
        return self.pipe[1]
    def run(self):
        for line in iter(self.pin.readline, ''):
            self.logger.log(self.level, line.rstrip())
        self.pin.close()
    def close(self):
        os.close(self.pipe[1])


def run(cmd, shell=False):
    logger.debug("executing: %s" % cmd)
    sublogger = loggers.get(cmd[0])
    errpipe = LoggerPipe(sublogger, logging.WARNING)
    outpipe = LoggerPipe(sublogger, logging.INFO)
    p = subprocess.Popen(cmd, stdout = outpipe, stderr = errpipe)
    p.wait()
    errpipe.close()
    outpipe.close()
    # output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=shell)
    # logger.debug("========command output========\n%s" % output.decode().rstrip())
    # logger.debug("==========end output==========")
    
