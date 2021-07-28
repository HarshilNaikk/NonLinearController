import matplotlib.pyplot as plt
import numpy as np

x, y, z, xdot, ydot, zdot, theta, phi, psi, thetadot, psidot, phidot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 

m = 2.0
Ix, Iy = 1.2416
Iz = 2.4832
d = 0.2
c = 0.01
g = 9.81
kt = np.diag([0.01, 0.01, 0.01])
kr = np.diag([0.001, 0.001, 0.001])
A1,A2,A3,A4,A5,A6 = np.diag([3,3])
A7 = np.diag([3,3,3,3])
dt = 0.1


x1 = np.array([x,y])
x2 = np.array([xdot, ydot])
x3 = np.array([theta, phi])
x4 = np.array([thetadot, phidot])
x5 = np.array([psi, alt])
x6 = np.array([psidot, altdot])
x7 = np.array([F1, F2, F3, F4])





# # Data for plotting
# t = np.arange(0.0, 2.0, 0.01)
# s = 1 + np.sin(2 * np.pi * t)

# fig, ax = plt.subplots()
# ax.plot(t, s)

# ax.set(xlabel='time (s)', ylabel='voltage (mV)',
#        title='About as simple as it gets, folks')
# ax.grid()

# fig.savefig("test.png")
# plt.show()