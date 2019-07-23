#!/usr/bin/env python3
import json
import os
import sys
import time
import argparse

caseroot = os.path.dirname(os.path.abspath(__file__))
case_vars = json.load(open(os.path.join(caseroot, "case.json")))
toolroot = case_vars["TOOLROOT"]
esmdroot = case_vars["ESMDROOT"]
LID=time.strftime("%y%m%d-%H%M%S")

sys.path.append(toolroot)
import loggers
import makedep
import shell
import shutil
import vc
parser = argparse.ArgumentParser(description='Create a new case with specified name')
parser.add_argument("action", default="build", action="store", help="type of action", choices=["clean", "test", "build"], type=str, nargs='?')
args = parser.parse_args(sys.argv[1:])


logger = loggers.get("build")
try:
    blddir = os.path.join(caseroot, 'bld')
    histdir = os.path.join(blddir, "history")
    objdir = os.path.join(blddir, "obj")
    gitsrcdir = os.path.join(blddir, "src")
    dirtygit = os.path.join(gitsrcdir, ".git")
    for d in [blddir, histdir, objdir, gitsrcdir]:
        if not os.path.exists(d):
            logger.info("%s does not exist, creating..." % d)
            os.makedirs(d)

    if not os.path.exists(dirtygit):
        logger.info("dirty git does not exist creating...")
        vc.git_init(gitsrcdir)
    bldlog = os.path.join(histdir, "bldlog.%s" % LID)
    loggers.add_file(bldlog)

    logger.info("copying files to temporary build directory")
    copied_src = set([])
    srcpath = [os.path.join(esmdroot, 'src'), os.path.join(caseroot, 'SourceMods')]
    platformsrc_path = []
    for srcdir in srcpath:
        subdir = os.path.join(srcdir, 'mods-%s' % case_vars["PLATFORM"])
        if os.path.exists(subdir):
            platformsrc_path.append(subdir)
    srcpath.extend(platformsrc_path)
    for srcdir in srcpath:
        for srcfile in os.listdir(srcdir):
            srcext = os.path.splitext(srcfile)[1]
            if srcext in [".c", ".h", ".i", ".s"]:
                shutil.copy2(os.path.join(srcdir, srcfile), gitsrcdir)
                copied_src.add(srcfile)
    logger.info("cleaning temporary build directory")
    os.chdir(gitsrcdir)
    for srcfile in os.listdir():
        srcext = os.path.splitext(srcfile)[1]
        if srcext in [".c", ".h", ".i", ".s"]:
            if srcfile  not in copied_src:
                os.remove(srcfile)
                logger.debug("removing %s", srcfile)

    logger.info("updating the dirty git")
    vc.git_add(add_all=True)
    vc.git_commit("auto commit for build %s" % LID)
    logger.info("generating dependencies...")
    
    incmap = makedep.gen_incmap(srcpath)
    deps = makedep.export_deps(incmap)
    
    deppath = os.path.join(objdir, "Depends")
    depfile = open(deppath, "w")
    depfile.write("%s\n" % deps)
    depfile.close()
    
    logger.info("dependencies written at %s", deppath)
    
    #VPATH=../../src:../../SourceMods/ TOOLROOT=../../../../tools/ USER_CPPDEFS="-DTEST" PLATFORM=gcc_arch
    logger.info("setting environment vars...")
    for var in ["USER_CPPDEFS", "PLATFORM"]:
        os.environ[var] = case_vars[var]
    os.environ["VPATH"] = ":".join(srcpath)
    os.environ["SRCPATH"] = ":".join(srcpath)
    os.environ["BLDDIR"] = os.path.join(caseroot, "bld")
    os.environ["HISTDIR"] = histdir
    os.environ["LID"] = LID
    os.environ["CASEROOT"] = caseroot
    os.chdir(objdir)
    ret = 0
    if args.action == "build":
        ret = shell.run(["make", "-f", os.path.join(toolroot, "Makefile")]).returncode
    elif args.action == "clean":
        ret = shell.run(["make", "-f", os.path.join(toolroot, "Makefile"), "clean"]).returncode
    if ret != 0:
        logger.error("build failed, exitting...")
        sys.exit(1)
    else:
        logger.info("done")
        if args.action == "build":
            shutil.copy2(os.path.join(objdir, "esmd"), os.path.join(histdir, "esmd-%s" % LID))
            shutil.copy2(os.path.join(objdir, "esmd"), os.path.join(blddir, "esmd"))
except Exception as e:
    logger.exception(e)
