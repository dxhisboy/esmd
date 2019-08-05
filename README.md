# ESMD

#### Description
Easy Scalable Molecular Dynamics Simulator.

#### Instructions

1. Use `create_case.py -P shenwei <casename>` to create a new case for Shenwei architecture.
2. Change directory to `cases/<casename>` and execute build.py (requires python3) to build esmd
3. Use `bsub -share_size 6144 -host_stack 16 -priv_size 4 -debug -I -n <#procs> -cgsp 64 -b -q <queue_name> -o reset.txt ./bld/esmd <nx> <ny> <nz> <cutoff> <nsteps>` to submit the program.

