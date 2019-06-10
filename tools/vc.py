import os
import loggers
import shell
logger = loggers.get('vc')
def init_git_args(worktree = None, gitdir = None, workdir = None):
    git_args = ["git"]
    if gitdir:
        git_args.append("--git-dir")
        git_args.append(gitdir)
    if worktree:
        git_args.append("--work-tree")
        git_args.append(worktree)
    if workdir:
        git_args.append("-C")
        git_args.append(workdir)
    return git_args
        
def git_init(path='.', **kwargs):
    git_args = init_git_args(**kwargs)
    git_args.append("init")
    git_args.append(path)
    logger.info('Initializing git repo at %s...', path)
    shell.run(git_args)

def git_add(path='.', add_all=False, **kwargs):
    git_args = init_git_args(**kwargs)
    git_args.append("add")
    git_args.append(path)
    shell.run(git_args)
    
def git_commit(message, **kwargs):
    git_args = init_git_args(**kwargs)
    git_args.append('commit')
    git_args.append('-m')
    git_args.append(message)
    shell.run(git_args)

def git_commit_hash(abbrev=True, **kwargs):
    git_args = init_git_args(**kwargs)
    git_args.append('show')
    git_args.append('-s')
    if abbrev:
        git_args.append('--format=%h')
    else:
        git_args.append('--format=%H')
    ret = shell.run(git_args, capture=True, quiet=True)
    if ret.returncode == 0:
        return ret.stdout.decode().strip()
    else:
        return None
