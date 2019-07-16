import sys
import re
import os
import subprocess
FUNC_RE = re.compile('(?P<deco>.*)\\s+(?P<name>[^\\s]*)\\s*\((?P<arg>.*)\)\\s*\x01')
COMML_RE = re.compile('//.*\n')
COMMR_RE = re.compile('/\*.*\*/')
DIR_RE = re.compile('#.*\n')
import array
#s = subprocess.check_output(["cpp", "-fpreprocessed", sys.argv[1]]).decode()
s = open(sys.argv[1]).read()
s = COMML_RE.sub("", s)
s = COMMR_RE.sub("", s)
s = s.replace("\\\n", "")
s = DIR_RE.sub("", s)
#print(s)
#print(s)
#s = open(sys.argv[1]).read()
buf = array.array('b')
nest = 0
indir = 0
for c in s:
    if c == '#':
        indir = 1
    if indir and c == '\n':
        indir = 0
        buf.append(0)
    if indir:
        continue
    if c == '{':
        nest += 1
        if nest == 1:
            buf.append(1)
        continue
    elif c == '}':
        nest -= 1
        if nest == 0:
            buf.append(0)
        continue

    if nest == 0 and indir == 0:
        if c == ';':
            buf.append(0)
        elif c in ['\n', '\t']:
            continue
        else:
            if nest == 0:
                buf.append(ord(c))
print(buf.tostring())
signatures = []
sp = buf.tostring().decode().split('\0')
for s in sp:
    m = FUNC_RE.match(s)
    if m:
        g = m.groupdict()
        if g['deco'].find('static') == -1 and g['deco'].find('inline') == -1:
            signatures.append("%s %s(%s);" % (g['deco'], g['name'], g['arg']))

header = os.path.splitext(sys.argv[1])[0] + ".h"
guard = os.path.splitext(sys.argv[1])[0].upper()
SIGSTART="//function signatures\n"
SIGEND="//end function signatures\n"
if not os.path.exists(header):
    hfile = open(header, "w")
    hfile.write("#ifndef %s_H_\n" % guard)
    hfile.write("#define %s_H_\n" % guard)
    hfile.write("#include <data.h>\n")
    hfile.write(SIGSTART)
    hfile.write(SIGEND)
    hfile.write("#endif\n")
    hfile.close()
hfile = open(header, "r")
hlines = hfile.readlines()
hfile.close()
hfile = open(header, "w")
insig = 0
for line in hlines:
    if line == SIGEND:
        insig = 0
    if insig == 0:
        hfile.write(line)
    if line == SIGSTART:
        insig = 1
        for sig in signatures:
            hfile.write(sig + '\n')
hfile.close()
