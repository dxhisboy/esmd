import argparse
import os
import sys
import shutil
import json

esmdroot = os.path.dirname(os.path.abspath(__file__))
toolroot = os.path.join(esmdroot, "tools")
sys.path.append(toolroot)

import loggers
import shell
import vc

parser = argparse.ArgumentParser(description='Create a new case with specified name')
parser.add_argument("casename", action="store", help="Name of case", type=str)
parser.add_argument("-p", "--proto", action="store", default=os.path.join(toolroot, "case.proto"), 
                    help="Create case as a copy of prototype case")
parser.add_argument("-P", "--platform", action="store", default=None,
                    help="Platform for the case")
parser.add_argument("-D", "--user-cppdefs", action="store", default="",
                    help="User defined preprocessing variables")
parser.add_argument("-R", "--keep-root", action="store", default="",
                    help="Keep root the same as prototype case")
args = parser.parse_args(sys.argv[1:])
args.proto = os.path.abspath(args.proto)
casesdir = os.path.join(esmdroot, "cases")
caseroot = os.path.join(esmdroot, "cases", args.casename)

logger = loggers.get("create_case")
try:
    if not os.path.exists(casesdir):
        logger.info("%s does not exist, creating..." % casesdir)
        
    logger.info("prototype case is: %s" % args.proto)
    logger.info("copying prototype to case: %s" % caseroot)
    shutil.copytree(args.proto, caseroot, ignore=shutil.ignore_patterns('bld', 'run', 'src', '.git'))
    casegitdir = os.path.join(caseroot, "SourceMods", ".git")
    protogitdir = os.path.join(args.proto, "SourceMods", ".git")
    if not os.path.exists(protogitdir):
        logger.info("prototype case has no SourceMods git, initializing...")
        #shell.run(["git", "init", os.path.join(caseroot, "SourceMods")], shell=False)
        vc.git_init(os.path.join(caseroot, "SourceMods"))
        vc.git_init(os.path.join(caseroot, "SourceMods"), ".dirty-git")
    else:
        logger.info("copying SourceMods git to %s" % casegitdir)
        shutil.copytree(protogitdir, casegitdir, copy_function=shutil.copy)
    logger.info("linking src to case")
    os.symlink(os.path.join(esmdroot, "src"), os.path.join(os.path.join(caseroot, "src")), True)

    logger.info("considering case.json")
    caseconf = os.path.join(os.path.join(caseroot, "case.json"))
    case_vars = {'PLATFORM' : 'gcc_arch', "USER_CPPDEFS" : ""}
    if os.path.exists(caseconf):
        logger.info("case.json exists, loading")
        case_vars.update(json.load(open(caseconf)))
    if args.platform:
        case_vars['PLATFORM'] = 'gcc_arch'
    if args.user_cppdefs:
        case_vars['USER_CPPDEFS'] += args.user_cppdefs
    if not args.keep_root:
        case_vars['ESMDROOT'] = esmdroot
        case_vars['TOOLROOT'] = toolroot
    logger.info("writing case.json...")
    output = open(caseconf, "w")
    json.dump(case_vars, output, indent=2)
    output.write("\n")
except Exception as e:
    logger.exception(e)
