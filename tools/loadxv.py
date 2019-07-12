import array
import numpy
# xv = array.array('d')
# bin = open("xv.bin", "rb").read()
# xv.frombytes(bin)

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
                                xi = (i) *  0.5
                                yi = (j) *  0.5
                                zi = (k) *  0.5
                                n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1
                                x.append(xi)
                                x.append(yi)
                                x.append(zi)
                                x.append(n)
    return numpy.asarray(x).reshape((131072, 4))
def create_fcc_by_offset(dens, nx, ny, nz):
    x = array.array('d')
    scale = (4 / dens) ** (1 / 3.)
    offset = [(0, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5)]
    for kk in range(0, nz):
        for jj in range(0, ny):
            for ii in range(0, nx):
                for o in offset:
                    xi = (ii + o[0]) * scale * 0.5
                    yi = (jj + o[1]) * scale * 0.5
                    zi = (kk + o[2]) * scale * 0.5
                    i = (ii + o[0]) * 2
                    j = (ii + o[1]) * 2
                    k = (ii + o[2]) * 2
                    
                    #n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1
                    x.append(xi)
                    x.append(yi)
                    x.append(zi)
                    x.append(0)
    return numpy.asarray(x).reshape((131072, 4))
xfo = create_fcc_by_offset(0.8442, 32, 32, 32)
xf = create_fcc(0.8442, 32, 32, 32)

