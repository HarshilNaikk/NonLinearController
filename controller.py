import matplotlib.pyplot as plt
import numpy as np
from numpy import sin as S 
from numpy import cos as C 
from numpy import tan as T 

x, y, z, xdot, ydot, zdot, theta, phi, psi, thetadot, psidot, phidot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 

m = 2.0
Ix, Iy = 1.2416, 1.2416
Iz = 2.4832
It = np.diag(Ix,Iy,Iz)
eta = np.array([phi,theta,psi])
etadot = np.array([phidot,thetadot,psidot])

d = 0.2
c = 0.01
g = 9.81
kt = np.diag([0.01, 0.01, 0.01])
kr = np.diag([0.001, 0.001, 0.001])
A1,A2,A3,A4,A5,A6 = np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),
A7 = np.diag([3,3,3,3])
dt = 0.1

g0 = ((F1+F2+F3+F4)/m)*np.array([[S(psi), C(psi)],[-C(psi), S(psi)]])
g1 = np.array([[1/Ix , S(phi)*T(theta)/Iy],[0, C(phi)/Iy]])
g2 = np.array([[C(phi)/(Iz*C(theta)) , 0],[0, C(phi)*C(theta)/m]])

phivec0 = np.array([S(phi)],[C(phi)*S(theta)])
phivec1 = np.array([[d*(F2-F4)],[d*(F3-F1)]])
phivec2 = np.array([[c*(F1-F2+F3-F4)],[F1+F2+F3+F4]])

Rt = np.array([C(phi)*C(psi), (S(phi)*S(theta)*C(psi) - C(phi)*S(psi)), (C(phi)*S(theta)*C(psi) + S(phi)*S(psi))], [C(theta)*S(psi), (S(phi)*S(theta)*S(psi) + C(phi)*C(psi)), (C(phi)*S(theta)*S(psi) - S(phi)*C(psi))], [-S(phi), S(phi)*C(theta), C(phi)*C(theta)])
Rr = np.array([[1,0,-S(theta)], [0,C(phi),C(theta)*S(phi)], [0,-S(phi),C(phi)*C(theta)]])
dRrbydphi = np.array([[0,0,0], [0,-S(phi)*phidot,C(theta)*C(phi)*phidot], [0,-C(phi)*phidot,-S(phi)*C(theta)*phidot]])
dRrbydphibydtheta = np.array([[0,0,-C(theta)*thetadot], [0,0,-S(theta)*S(phi)*thetadot], [0,0,-C(phi)*S(theta)*thetadot]])

fx, fy, fz = 0.0, 0.0, -g
fphi,ftheta,fpsi = 0.0, 0.0, 0.0
farr1 = np.array([[fx],[fy],[fz]])
farr2 = np.array([[fphi],[ftheta],[fpsi]])
farr2 = -np.linalg.inv(np.dot(It,Rr))*(It.dot(phidot*dRrbydphi+thetadot*dRrbydphibydtheta).dot(etadot) - np.cross(Rr.dot(etadot), It.dot(Rr).dot(etadot))) + np.array([[(c/Iz)*C(phi)*T(theta)*(F1-F2+F3-F4)],[(-c/Iz)*S(phi)*(F1-F2+F3-F4)],[(d/Iy)*(S(phi)/T(theta))*(F3-F1)]])

f0 = np.array([[fx],[fy]])
f1 = np.array([[fphi],[fpsi]])
f2 = np.array([[fpsi],[fz]])


x1 = np.array([x,y])
x2 = np.array([xdot, ydot])
x3 = np.array([theta, phi])
x4 = np.array([thetadot, phidot])
x5 = np.array([psi, alt])
x6 = np.array([psidot, altdot])
x7 = np.array([F1, F2, F3, F4])


def integrate(self,q1,q2):
    return q2*self.dt + q1

self.x1dot = self.x2
self.x2dot = self.f0 + self.g0*self.phivec0
self.x3dot = self.x4
self.x4dot = self.f1 + self.g1*self.phivec1
self.x5dot = self.x6
self.x6dot = self.f2 + self.g2*self.phivec2
self.x7dot = self.u

self.x7 = integrate(self.x7,self.x7dot)
self.x6 = integrate(self.x6,self.x6dot)
self.x5 = integrate(self.x5,self.x5dot)
self.x4 = integrate(self.x4,self.x4dot)
self.x3 = integrate(self.x3,self.x3dot)
self.x2 = integrate(self.x2,self.x2dot)
self.x1 = integrate(self.x1,self.x1dot)




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