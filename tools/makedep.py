import re
import os
import json
import time
INCLUDE_RE = re.compile('\\s*#include\\s*((<(?P<sysinc>.*)>)|("(?P<usrinc>.*)"))')
TEMPLATE_RE = re.compile('\\s*#define\\s*TEMPLATE\\s*((<(?P<sysinc>.*)>)|("(?P<usrinc>.*)"))')
src_suffix = set([".c", ".h"])

def target_name(f):
    if f.endswith(".c"):
        base = os.path.basename(f)
        return os.path.splitext(base)[0] + ".o"
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
    for line in open(find_file(f, srcpath)):
        m = INCLUDE_RE.match(line) or TEMPLATE_RE.match(line)
        if m:
            groups = m.groupdict()
            inc = groups['sysinc'] or groups['usrinc']
            inc_full = find_file(inc, srcpath)

            if inc_full:
                incmap[f].add(inc)
                find_includes(inc, incmap, srcpath)
                incmap[f].update(incmap[inc])

def gen_incmap(srcpath, initial={}):
    incmap = dict(initial.items())

    for path in srcpath:
        for f in os.listdir(path):
            f_full = os.path.join(path, f)
            find_includes(f_full, incmap, srcpath)
    return incmap

def export_deps(incmap, target_name = target_name, objvar = "OBJS"):
    deps = []
    objs = []
    for k, v in incmap.items():
        obj = target_name(k)
        if obj is not None:
            deps.append("%s: %s" % (obj, " ".join(v)))
            objs.append(obj)
    objline = ""
    if objvar:
        objline = "%s=%s\n" % (objvar, " ".join(objs))
    return objline + os.linesep.join(deps)

def save_incmap(incmap, path):
    incmap_encodable = dict(map(lambda k: (k[0], list(k[1])), incmap.items()))
    record = {"stamp": time.time(), "map": incmap_encodable}
    #encode = json.JSONEncoder().encode(record)
    #open(path, 'w').write(str(encode))
    json.dump(record, open(path, 'w'))
def check_existing_incmap(path, srcpath):
    if not os.path.exists(path):
        return {}
    encode = open(path, 'r').read()
    record = json.JSONDecoder().decode(encode)
    stamp = record["stamp"]
    incmap_orig = record["map"]
    incmap_new = {}
    for k in incmap_orig:
        if find_file(k, srcpath):
            if os.stat(k).st_mtime < stamp:
                incmap_new[k] = list(incmap_orig[k])
    return incmap_new
