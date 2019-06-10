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
    def join(self):
        super().join()
    def close(self):
        os.close(self.pipe[1])
def run(cmd, shell=False, capture=False, quiet=False, text=True):
    if not quiet:
        logger.debug("executing: %s" % cmd)
    p = None
    if not capture:
        if not quiet:
            sublogger = loggers.get(cmd[0])
            errpipe = LoggerPipe(sublogger, logging.WARNING)
            outpipe = LoggerPipe(sublogger, logging.DEBUG)
        else:
            errpipe = subprocess.DEVNULL
            outpipe = subprocess.DEVNULL
        
        p = subprocess.run(cmd, stdout = outpipe, stderr = errpipe, text=text, universal_newlines=True)
        if not quiet:
            errpipe.close()
            outpipe.close()
            errpipe.join()
            outpipe.join()
    else:
        p = subprocess.run(cmd, capture_output=True)
    return p
    # output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=shell)
    # logger.debug("========command output========\n%s" % output.decode().rstrip())
    # logger.debug("==========end output==========")
     
