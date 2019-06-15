import subprocess
import sys
import os
def add_args(parser):
    pass
def dispatch(cmd, args=None, logfile=None, **kwargs):
    if logfile:
        tee = subprocess.Popen(["tee", logfile], stdout=sys.stdout, stderr=sys.stdout, stdin=subprocess.PIPE)
        proc = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=tee.stdin)
        proc.wait()
        tee.stdin.close()
        tee.wait()
    else:
        proc = subprocess.Popen(cmd, stderr=sys.stderr, stdout=sys.stdout)
        proc.wait()
