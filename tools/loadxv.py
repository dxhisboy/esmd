import array
import numpy
xv = array.array('d')
bin = open("xv.bin", "rb").read()
xv.frombytes(bin)

# x = array.array('d')
# v = array.array('d')
# for i in range(0, len(xv), 6):
#     x.append(xv[i])
#     x.append(xv[i + 1])
#     x.append(xv[i + 2])
#     v.append(xv[i + 3])
#     v.append(xv[i + 4])
#     v.append(xv[i + 5])

# x = numpy.asarray(x).reshape((131072, 3))
# v = numpy.asarray(v).reshape((131072, 3))
def create_fcc(dens, nx, ny, nz):
    x = array.array('d')
    scale = (4 / dens) ** (1 / 3.)
    offset = [(0, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5)]
    for kk in range(0, nz * 2, 8):
        for jj in range(0, ny * 2, 8):
            for ii in range(0, nx * 2, 8):
                for k in range(kk, kk + 8):
                    for j in range(jj, jj + 8):
                        for i in range(ii, ii + 8):
                            if ((i + j + k) % 2 == 0):
                                xi = (i) * scale * 0.5
                                yi = (j) * scale * 0.5
                                zi = (k) * scale * 0.5
                                n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1
                                x.append(xi)
                                x.append(yi)
                                x.append(zi)
                                x.append(n)
    return numpy.asarray(x).reshape((131072, 4))

IA=16807
IM=2147483647
AM=(1.0/IM)
IQ=127773
IR=2836
MASK=123459876

def random(l):
    seed = l[0]
    k = int(seed / IQ)
    seed = IA * (seed - k * IQ) - IR * k
    if seed < 0:
        seed += IM
    l[0] = seed
    return AM * seed
def create_fcc_by_offset(dens, nx, ny, nz):
    x = array.array('d')
    scale = (4 / dens) ** (1 / 3.)
    offset = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    for kk in range(0, nz):
        for jj in range(0, ny):
            for ii in range(0, nx):
                for o in offset:
                    xi = (ii + o[0]) * scale
                    yi = (jj + o[1]) * scale
                    zi = (kk + o[2]) * scale
                    i = xi / scale * 2
                    j = yi / scale * 2
                    k = zi / scale * 2
                    
                    n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1
                    l = [n]
                    
                    x.append(xi)
                    x.append(yi)
                    x.append(zi)
                    for m in range(5):
                        random(l)
                    vx = random(l)
                    for m in range(5):
                        random(l)
                    vy = random(l)
                    for m in range(5):
                        random(l)
                    vz = random(l)
                    x.append(vx)
                    x.append(vy)
                    x.append(vz)
    return numpy.asarray(x).reshape((131072, 4))
xo = create_fcc_by_offset(0.8442, 32, 32, 32)
xf = create_fcc(0.8442, 32, 32, 32)

xo = numpy.asarray(sorted(xo.tolist(), key=lambda x: x[0] * 10000 + x[1] * 100 + x[2]))
xf = numpy.asarray(sorted(xf.tolist(), key=lambda x: x[0] * 10000 + x[1] * 100 + x[2]))

