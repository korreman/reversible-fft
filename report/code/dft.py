# naive reference implementation of DFT for comparison with janus FFT
import math
from matplotlib import pyplot as plt

def rotate(x, angle):
    s = math.sin(angle)
    c = math.cos(angle)
    return (c * x[0] - s * x[1], s * x[0] + c * x[1])

def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def dft(signal):
    l = len(signal)
    result = [(0,0) for x in signal]

    for k in range(l):
        for n, x in zip(range(len(signal)), signal):
            result[k] = add(result[k], rotate(x, k * n / len(signal) * math.pi * 2))
    return result

signal = [0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080,0,25080,46341,60547,65536,60547,46341,25080,0,-25080,-46341,-60547,-65536,-60547,-46341,-25080]
signal = [(s, 0) for s in signal]

result = dft(signal)
yr = [r for r, i in result]
yi = [i for r, i in result]

x = range(256)

y = [math.sqrt(r**2 + i**2) for r, i in zip(yr, yi)]
plt.plot(x, y)
plt.show()

ax = plt.axes(projection = '3d')
ax.set_proj_type('ortho')

ax.plot3D(x, yr, yi)
plt.show()
