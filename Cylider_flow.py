import numpy as np
import matplotlib.pyplot as plt

# Problem Specifications:
error, L, H, B, La, Re = 1e-03, 18.0, 10.0, 1.0, 7.5, 24

# Domain:
del_t, alpha, iim, jim = 0.02, 0.2, 180, 84
ire, jre = iim - 1, jim - 1
dx, dy, u_inf = (L / B) / ire, (H / B) / jre, 1
Rex, Rey = dx / (Re * dy), dy / (Re * dx)
als, ale, bls, ble = round(((ire - 1) / L) * (La - (B / 2))), round(((ire - 1) / L) * (La + (B / 2))), round(((jre - 1) / H) * ((H / 2) - (B / 2))), round(((jre - 1) / H) * ((H / 2) + (B / 2)))

def matrix_initialize(zmax):
    b = np.empty([zmax])
    d = np.empty([zmax])
    a = np.empty([zmax])
    c = np.empty([zmax])
    fv = np.empty([zmax])
    return (b, d, a, c, fv)

def final(b, d, a, c, zmax, fv):
    for z in range(1, zmax):
        d[z] -= b[z] * (a[z - 1] / d[z - 1])
        c[z] -= b[z] * (c[z - 1] / d[z - 1])
        fv[zmax - 1] = c[zmax - 1] / d[zmax - 1]
    for z in range(zmax - 2, -1, -1):
        fv[z] = (c[z] - a[z] * fv[z + 1]) / d[z]
    return fv

def bottom():
    for i in range(1, ire):
        u[0, i] = -u[1, i]
        v[0, i] = 0.0

def top():
    for i in range(1, ire):
        u[jim - 1, i] = -u[jre - 1, i]
        v[jre - 1, i] = 0.0

def left():
    for j in range(0, jim):
        u[j, 0] = u_inf
        v[j, 0] = -v[j, 1]

def right():
    for j in range(0, jim):
        u[j, ire - 1] = u[j, ire - 2]
        v[j, ire] = v[j, ire - 1]
    dp[j, ire] = 0

def cylinder():
    for i in range(als, ale - 1):
        u[bls, i] = -u[bls - 1, i]
        v[bls - 1, i] = 0.0
    for i in range(als, ale - 1):
        u[ble - 1, i] = -u[ble, i]
        v[ble - 1, i] = 0.0
    for j in range(bls - 1, ble):
        u[j, als - 1] = 0.0
    for j in range(bls - 1, ble):
        v[j, als] = -v[j, als - 1]
    for j in range(bls - 1, ble):
        u[j, ale - 1] = 0.0
    for j in range(bls - 1, ble):
        v[j, ale - 1] = -v[j, ale]

aeu = lambda j, i: (u[j, i] + u[j, i + 1]) * dy / 4 - (Rey)
awu = lambda j, i: -(u[j, i - 1] + u[j, i]) * dy / 4 - (Rey)
anu = lambda j, i: (v[j, i] + v[j, i + 1]) * dx / 4 - (Rex)
asu = lambda j, i: -(v[j - 1, i] + v[j - 1, i + 1]) * dx / 4 - (Rex)
bu = lambda j, i: (-(dx * dy) / del_t) * u[j, i]
apu = lambda j, i: (p[j, i + 1] - p[j, i]) * dy
au = lambda j, i: 0.25 * (u[j, i] + u[j, i + 1]) * dy - 0.25 * (u[j, i - 1] + u[j, i]) * dy + 0.25 * (v[j, i] + v[j, i + 1]) * dx - 0.25 * (v[j - 1, i] + v[j - 1, i + 1]) * dx + 2 * (Rey) + 2 * (Re)
aeu = lambda j, i: (u[j, i] + u[j, i + 1]) * dy / 4 - (Rey)
awu = lambda j, i: -(u[j, i - 1] + u[j, i]) * dy / 4 - (Rey)
anu = lambda j, i: (v[j, i] + v[j, i + 1]) * dx / 4 - (Rex)
asu = lambda j, i: -(v[j - 1, i] + v[j - 1, i + 1]) * dx / 4 - (Rex)
bu = lambda j, i: (-(dx * dy) / del_t) * u[j, i]
apu = lambda j, i: (p[j, i + 1] - p[j, i]) * dy
au = lambda j, i: 0.25 * (u[j, i] + u[j, i + 1]) * dy - 0.25 * (u[j, i - 1] + u[j, i]) * dy + 0.25 * (v[j, i] + v[j, i + 1]) * dx - 0.25 * (v[j - 1, i] + v[j - 1, i + 1]) * dx + 2 * (Rey) + 2 * (Re)

aev = lambda j, i: (u[j, i] + u[j + 1, i]) * dy / 4 - (Rey)
awv = lambda j, i: -(u[j, i - 1] + u[j + 1, i - 1]) * dy / 4 - (Rey)
anv = lambda j, i: (v[j, i] + v[j + 1, i]) * dx / 4 - (Rex)
asv = lambda j, i: -(v[j, i] + v[j - 1, i]) * dx / 4 - (Rex)
bv = lambda j, i: (-(dx * dy) / del_t) * v[j, i]
apv = lambda j, i: (p[j + 1, i] - p[j, i]) * dx
av = lambda j, i: 0.25 * (u[j, i] + u[j + 1, i]) * dy - 0.25 * (u[j, i - 1] + u[j + 1, i - 1]) * dy + 0.25 * (v[j, i] + v[j + 1, i]) * dx - 0.25 * (v[j, i] + v[j - 1, i]) * dx + 2 * (Rey) + 2 * (Re)

def u_correct():
    b, d, a, c, fv = matrix_initialize(iim)
    for j in range(1, jre):
        cylinder()
        for i in range(1, ire):
            b[0], d[0], a[0], c[0] = 0, 1, 0, u_inf
            b[-1], d[-1], a[-1], c[-1] = -1, 1, 0, 0
            b[i], d[i], a[i], c[i] = awu(j, i), au(j, i) + ((dx * dy) / del_t), aeu(j, i), -bu(j, i) - apu(j, i) - asu(j, i) * u[j - 1, i] - anu(j, i) * u[j + 1, i]
            fv[i] = u[j, i]
        u[j, :] = final(b, d, a, c, iim, fv)
    bottom()
    top()
    b, d, a, c, fv = matrix_initialize(jim)
    for i in range(1, ire):
        for j in range(1, jre):
            b[0], d[0], a[0], c[0] = 0, 1, 0, u_inf
            b[-1], d[-1], a[-1], c[-1] = -1, 1, 0, 0
            b[i], d[i], a[i], c[i] = awu(j, i), av[j, i] + ((dx * dy) / del_t), aev(j, i), -bv(j, i) - apv(j, i) - asv(j, i) * v[j - 1, i] - anv(j, i) * v[j + 1, i]
            fv[i] = v[j, i]
        v[j, :] = final(b, d, a, c, jim, fv)
    left()
    right()

def p_tridiag():
    global p
    b, d, a, c, fv = matrix_initialize(iim)
    for j in range(1, jre):
        for i in range(1, ire):
            b[i], d[i], a[i], c[i] = awu(j, i) - aeu(j, i) - avu(j, i) - anu(j, i), 2 * (dx * dy) / del_t, aev(j, i) - avu(j, i) - anu(j, i), -aev(j, i) - avu(j, i) - anu(j, i)
            fv[i] = (u[j + 1, i] - u[j, i]) * dx
        p[j, :] = final(b, d, a, c, iim, fv)
    bottom()
    top()

def dp_solution():
    global p
    global dp
    for i in range(1, ire - 1):
        for j in range(1, jre - 1):
            dp[j, i] = (dx / (((dx * dy) / del_t) + av[j, i])) * (dp[j, i] - dp[j + 1, i])

def algorithm():
    global a
    global dp
    global p
    global u
    global v
    global error
    global z
    global b_p
    while a > error:
    left()
    right()
    top()
    bottom()
    cylinder()
    u_correct()
    p_tridiag()
    dp_solution()
    print(z)
    z = z + 1
    a = np.amax(b_p)

u[bls-1:ble+1,als-1:ale+1]=0
v[bls-1:ble+1,als-1:ale+1]=0

# Plotting
x = np.linspace(dx, (L / B) - dx, ire - 1)
y = np.linspace(dy, (H / B) - dy, jre - 1)
X, Y = np.meshgrid(x, y)
plt.streamplot(X, Y, u[1:jre, 1:ire], v[1:jre, 1:ire], density=2, color='k')
plt.title("streamline")
plt.xlabel("x")
plt.ylabel("y")

# Plotting resultant velocity
skip = (slice(None, None, 4), slice(None, None, 4))
plt.quiver(X, Y, u[1:jre, 1:ire], v[1:jre, 1:ire], color='k')
plt.show()  
