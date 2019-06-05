import os
import loggers
import shell
logger = loggers.get('vc')
def git_init(path='.', gitdir = None):
    git_arg = ["git"]
    if gitdir:
        git_arg.append("--git-dir")
        git_arg.append(gitdir)
    git_arg.append("init")
    git_arg.append(path)
    logger.info('Initializing git repo at %s...', path)
    shell.run(git_arg)
