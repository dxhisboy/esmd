import re
import os
INCLUDE_RE = re.compile('\\s*#include\\s*((<(?P<sysinc>.*)>)|("(?P<usrinc>.*)"))')
src_suffix = set([".c", ".h"])

def target_name(f):
    if f.endswith(".c"):
        return os.path.basename(f) + ".o"
    return None

def find_file(f, srcpath):
    for path in srcpath:
        if os.path.exists(os.path.join(path, f)):
            return os.path.join(path, f)
    return None

def find_includes(f, incmap, srcpath):
    if os.path.splitext(f)[1] not in src_suffix:
        return
    if f in incmap:
        return
    incmap[f] = set([])
    for line in open(f):
        m = INCLUDE_RE.match(line)
        if m:
            groups = m.groupdict()
            inc = groups['sysinc'] or groups['usrinc']
            inc_full = find_file(inc, srcpath)

            if inc_full:
                incmap[f].add(inc_full)
                find_includes(inc_full, incmap, srcpath)
                incmap[f].update(incmap[inc_full])

def gen_incmap(srcpath):
    incmap = {}

    for path in srcpath:
        for f in os.listdir(path):
            f_full = os.path.join(path, f)
            find_includes(f_full, incmap, srcpath)
    return incmap

def export_dep(incmap, target_name = target_name):
    ret = []
    for k, v in incmap.items():
        obj = target_name(k)
        if obj is not None:
            ret.append("%s: %s" % (obj, " ".join(v)))
    return os.linesep.join(ret)

# print(export_dep(gen_incmap(["../src"])))
# gen_dep(["../src"])
# for f in os.listdir("../src"):
#     find_includes(os.path.join("../src", f), incmap, srcpath)
# print(incmap)
# for k, v in incmap.items():
#     obj = target_name(k)
#     if obj is not None:
#         print("%s: %s" % (obj, " ".join(v)))
