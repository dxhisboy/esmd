import json
import os
import sys
import time

caseroot = os.path.dirname(os.path.abspath(__file__))
case_vars = json.load(open(os.path.join(caseroot, "case.json")))
toolroot = case_vars["TOOLROOT"]
esmdroot = case_vars["ESMDROOT"]
LID=time.strftime("%y%m%d-%H%M%S")

sys.path.append(toolroot)
import loggers
import makedep
import shell

logger = loggers.get("build")

blddir = os.path.join(caseroot, 'bld')
histdir = os.path.join(blddir, "history")
objdir = os.path.join(blddir, "obj")
for d in [blddir, histdir, objdir]:
    if not os.path.exists(d):
        logger.info("%s does not exist, creating..." % d)
        os.makedirs(d)

bldlog = os.path.join(histdir, "bldlog.%s" % LID)
loggers.add_file(bldlog)
logger.info("generating dependencies...")

srcpath = [os.path.join(caseroot, 'SourceMods'), os.path.join(esmdroot, 'src')]
incmap = makedep.gen_incmap(srcpath)
deps = makedep.export_deps(incmap)

deppath = os.path.join(objdir, "Depends")
depfile = open(deppath, "w")
depfile.write("%s\n" % deps)
depfile.close()

logger.info("dependencies written at %s", deppath)

#VPATH=../../src:../../SourceMods/ TOOLROOT=../../../../tools/ USER_CPPDEFS="-DTEST" PLATFORM=gcc_arch
logger.info("setting environment vars...")
for var in ["TOOLROOT", "USER_CPPDEFS", "PLATFORM"]:
    os.environ[var] = case_vars[var]
os.environ["VPATH"] = ":".join(srcpath)
os.environ["BLDDIR"] = os.path.join(caseroot, "bld")
os.environ["HISTDIR"] = histdir
os.environ["LID"] = LID
os.chdir(objdir)
shell.run(["make", "-f", os.path.join(toolroot, "Makefile")])

