# naive reference implementation of DFT for comparison with janus FFT
import math
from matplotlib import pyplot as plt

def plot3D(arr, title):
    x = range(len(arr))
    yr = [r for r, i in arr]
    yi = [i for r, i in arr]
    ax = plt.axes(projection = '3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
    ax.set_ylim3d(-2e6, 2e6)
    ax.set_zlim3d(-2e6, 2e6)
    ax.plot3D(x, yr, yi)
    plt.title(title)
    plt.show()

def plot2D(arr):
    x = range(len(arr))
    yr = [r for r, i in arr]
    yi = [i for r, i in arr]
    y = [math.sqrt(r**2 + i**2) for r, i in zip(yr, yi)]
    plt.plot(x, y)
    plt.show()

twiddle_m = [0,-805,-1609,-2414,-3220,-4027,-4835,-5644,-6455,-7268,-8084,-8901,-9722,-10545,-11372,-12202,-13036,-13875,-14717,-15564,-16416,-17274,-18137,-19006,-19881,-20762,-21651,-22546,-23450,-24361,-25281,-26209,-27146,-28093,-29051,-30018,-30997,-31987,-32989,-34003,-35030,-36071,-37126,-38196,-39281,-40383,-41501,-42636,-43790,-44963,-46156,-47370,-48605,-49863,-51145,-52452,-53785,-55144,-56532,-57950,-59399,-60880,-62396,-63947]
twiddle_i = [0,1608,3215,4821,6423,8022,9616,11204,12785,14359,15923,17479,19024,20557,22078,23586,25079,26557,28020,29465,30893,32302,33692,35061,36409,37736,39039,40319,41575,42806,44011,45189,46340,47464,48558,49624,50660,51665,52639,53581,54491,55368,56212,57022,57797,58538,59243,59913,60547,61144,61705,62228,62714,63162,63571,63943,64276,64571,64826,65043,65220,65358,65457,65516]

def rotate_lifting(x, k, k_step):
    m = twiddle_m[k * k_step]
    i = twiddle_i[k * k_step]
    (xr, xi) = x
    xr += (m * xi) >> 16
    xi += (i * xr) >> 16
    xr += (m * xi) >> 16
    return (xr, xi)

def rotate_lifting_backward(x, k, k_step):
    m = twiddle_m[k * k_step]
    i = twiddle_i[k * k_step]
    (xr, xi) = x
    xr -= (m * xi) >> 16
    xi -= (i * xr) >> 16
    xr -= (m * xi) >> 16
    return (xr, xi)

def rotate(x, angle):
    s = math.sin(angle)
    c = math.cos(angle)
    return (c * x[0] - s * x[1], s * x[0] + c * x[1])

def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])

# NAIVE DFT
def dft(signal):
    l = len(signal)
    result = [(0,0) for x in signal]

    for k in range(l):
        for n, x in zip(range(len(signal)), signal):
            result[k] = add(result[k], rotate(x, k * n / len(signal) * math.pi * 2))
    return result

# NAIVE FFT
def fft_simple(signal):
    N = len(signal)
    halfN = N // 2
    if N == 1:
        return signal

    evens = [signal[2*n] for n in range(halfN)]
    odds = [signal[2*n + 1] for n in range(halfN)]

    evens = fft(evens)
    odds = fft(odds)

    for k in range(halfN):
        o = rotate(odds[k], 2 * math.pi * k/N)
        odds[k] = sub(evens[k], o)
        evens[k] = add(evens[k], o)

    return evens + odds

# LAYERED FFT
def scramble(signal):
    N = len(signal)
    halfN = N // 2
    if N == 1:
        return signal
    evens = [signal[2*n] for n in range(halfN)]
    odds = [signal[2*n + 1] for n in range(halfN)]
    evens = scramble(evens)
    odds = scramble(odds)
    return evens + odds

def fft_rotate(signal, N):
    halfN = N // 2
    quarterN = halfN // 2
    if halfN > 1:
        num_sections = len(signal) // N
        for section in range(num_sections):
            offset = halfN + section * N
            for k in range(1, quarterN): # set to halfN for naive
                odd_idx = k + offset
                dual_idx = halfN - k + offset
                signal[odd_idx] = rotate_lifting(signal[odd_idx], k, num_sections)
                # rotation without using lifting steps
                #signal[odd_idx] = rotate(signal[odd_idx], 2 * math.pi * k/N)
                #tmp = rotate(signal[dual_idx], -2 * math.pi * k/N)
                tmp = rotate_lifting_backward(signal[dual_idx], k, num_sections)
                signal[dual_idx] = (-tmp[0], -tmp[1])
            quarter_idx = offset + quarterN
            tmp = signal[quarter_idx]
            signal[quarter_idx] = (-tmp[1], tmp[0])
    return signal

def fft_convolve(signal, N):
    halfN = N // 2
    for section in range(len(signal) // N):
        for k in range(halfN):
            even_idx = k + section * N
            odd_idx = even_idx + halfN
            o = signal[odd_idx]
            signal[odd_idx] = sub(signal[even_idx], o)
            signal[even_idx] = add(signal[even_idx], o)
    return signal

def fft_layered(signal):
    signal = scramble(signal)
    #plot3D(signal, "Scrambled")
    N = 2
    while N <= len(signal):
        signal = fft_rotate(signal, N)
        #plot3D(signal, "Rotated")
        signal = fft_convolve(signal, N)
        #plot3D(signal, "Convolved")
        N *= 2
    return signal

def gen_sine(freq, len):
    return [math.floor(2**16 * math.sin(1/freq * 2 * math.pi * x)) for x in range(len)]

signal = [(x + y, 0) for x, y in zip(gen_sine(3, 256), gen_sine(7, 256))]
#print("signal:", [x for (x, y) in signal])
#plot3D(signal, "Initial")
result = fft_layered(signal)
plot2D(result)
