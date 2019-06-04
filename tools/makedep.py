import re
import os
import json
import time
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

def gen_incmap(srcpath, initial={}):
    incmap = dict(initial.items())

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
