import matplotlib.pyplot as plt
import numpy as np
from numpy import sin as S 
from numpy import cos as C 
from numpy import tan as T 

class controller:
    def __init__(self):
        self.x, self.y, self.z, self.xdot, self.ydot, self.zdot, self.theta, self.phi, self.psi, self.thetadot, self.psidot, self.phidot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 
        self.xd, self.yd, self.zd, self.psid = 0.0,0.0,0.0,

        self.m = 2.0
        self.Ix, self.Iy = 1.2416, 1.2416
        self.Iz = 2.4832
        self.It = np.diag(self.Ix,self.Iy,self.Iz)
        self.eta = np.array([self.phi,self.theta,self.psi])
        self.etadot = np.array([self.phidot,self.thetadot,self.psidot])
        self.F1, self.F2, self.F3, self.F4 = self.m*self.g/4,self.m*self.g/4,self.m*self.g/4,self.m*self.g/4,
        
        self.d = 0.2
        self.c = 0.01
        self.g = 9.81
        self.kt = np.diag([0.01, 0.01, 0.01])
        self.kr = np.diag([0.001, 0.001, 0.001])
        self.A1,self.A2,self.A3,self.A4,self.A5,self.A6 = np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),
        self.A7 = np.diag([3,3,3,3])
        self.dt = 0.1
        self.t = 0

        self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u = 0,0,0,0,0,0,0

    def update(self):
        self.g0 = ((self.F1+self.F2+self.F3+self.F4)/self.m)*np.array([[S(self.psi), C(self.psi)],[-C(self.psi), S(self.psi)]])
        self.g1 = np.array([[1/self.Ix , S(self.phi)*T(self.theta)/self.Iy],[0, C(self.phi)/self.Iy]])
        self.g2 = np.array([[C(self.phi)/(self.Iz*C(self.theta)) , 0],[0, C(self.phi)*C(self.theta)/self.m]])

        self.phivec0 = np.array([S(self.phi)],[C(self.phi)*S(self.theta)])
        self.phivec1 = np.array([[self.d*(self.F2-self.F4)],[self.d*(self.F3-self.F1)]])
        self.phivec2 = np.array([[self.c*(self.F1-self.F2+self.F3-self.F4)],[self.F1+self.F2+self.F3+self.F4]])

        self.Rt = np.array([C(self.phi)*C(self.psi), (S(self.phi)*S(self.theta)*C(self.psi) - C(self.phi)*S(self.psi)), (C(self.phi)*S(self.theta)*C(self.psi) + S(self.phi)*S(self.psi))], [C(self.theta)*S(self.psi), (S(self.phi)*S(self.theta)*S(self.psi) + C(self.phi)*C(self.psi)), (C(self.phi)*S(self.theta)*S(self.psi) - S(self.phi)*C(self.psi))], [-S(self.phi), S(self.phi)*C(self.theta), C(self.phi)*C(self.theta)])
        self.Rr = np.array([[1,0,-S(self.theta)], [0,C(self.phi),C(self.theta)*S(self.phi)], [0,-S(self.phi),C(self.phi)*C(self.theta)]])
        self.dRrbydphi = np.array([[0,0,0], [0,-S(self.phi)*self.phidot,C(self.theta)*C(self.phi)*self.phidot], [0,-C(self.phi)*self.phidot,-S(self.phi)*C(self.theta)*self.phidot]])
        self.dRrbydphibydtheta = np.array([[0,0,-C(self.theta)*self.thetadot], [0,0,-S(self.theta)*S(self.phi)*self.thetadot], [0,0,-C(self.phi)*S(self.theta)*self.thetadot]])

        self.fx, self.fy, self.fz = 0.0, 0.0, -self.g
        self.fphi,self.ftheta,self.fpsi = 0.0, 0.0, 0.0
        self.farr1 = np.array([[self.fx],[self.fy],[self.fz]])
        self.farr2 = np.array([[self.fphi],[self.ftheta],[self.fpsi]])
        self.farr2 = -np.linalg.inv(np.dot(self.It,self.Rr))*(self.It.dot(self.phidot*self.dRrbydphi+self.thetadot*self.dRrbydphibydtheta).dot(self.etadot) - np.cross(self.Rr.dot(self.etadot), self.It.dot(self.Rr).dot(self.etadot))) + np.array([[(self.c/self.Iz)*C(self.phi)*T(self.theta)*(self.F1-self.F2+self.F3-self.F4)],[(-self.c/self.Iz)*S(self.phi)*(self.F1-self.F2+self.F3-self.F4)],[(self.d/self.Iy)*(S(self.phi)/T(self.theta))*(self.F3-self.F1)]])

        self.f0 = np.array([[self.fx],[self.fy]])
        self.f1 = np.array([[self.fphi],[self.fpsi]])
        self.f2 = np.array([[self.fpsi],[self.fz]])


        self.x1 = np.array([self.x,self.y])
        self.x2 = np.array([self.xdot, self.ydot])
        self.x3 = np.array([self.theta, self.phi])
        self.x4 = np.array([self.thetadot, self.phidot])
        self.x5 = np.array([self.psi, self.z])
        self.x6 = np.array([self.psidot, self.zdot])
        self.x7 = np.array([self.F1, self.F2, self.F3, self.F4])


        self.x5d = np.array([0, 0])
        self.x5ddot = np.array([0, 0])
        self.x1ddot, self.x1d = np.array([1, C(self.t)]), np.array([self.t, S(self.t)])

        self.x1dot = self.x2
        self.x2dot = self.f0 + self.g0*self.phivec0
        self.x3dot = self.x4
        self.x4dot = self.f1 + self.g1*self.phivec1
        self.x5dot = self.x6
        self.x6dot = self.f2 + self.g2*self.phivec2
        self.x7dot = self.u

        self.jvec = np.array([J1,J2])
        #define Js, gblock,varrs, make sure about A* is dot product or not, define dot products in calc_input

    # def desiredVals(self,t):
    #    return np.array([1, np.cos(t)]), np.array([t, np.sin(t)])

    def calc_input(self):
        self.v1 = self.A1*(self.x1d-self.x1) + self.x1ddot
        v1dot = self.derivative(self.v1,self.v1o)
        self.v2 = np.linalg.inv(self.g0)*( (self.x1d - self.x1) + self.A2*(self.v1 - self.x2) + v1dot - self.f0 )
        self.v3 = np.linalg.inv(self.J0)*( np.transpose(self.g0)*(self.v1 - self.x2) + self.A3*(self.v2 - self.phivec0) + v2dot )
        self.v4 = np.linalg.inv(self.g1)*( np.transpose(self.J0)*(self.v2 - self.phivec0) + self.A4*(self.v3 - self.x4) + v3dot - self.f1 )
        self.v5 = self.A5*(self.x5d - self.x5) + self.x5ddot
        self.v6 = np.linalg.inv(self.g2)*( (self.x5d - self.x5) + self.A6*(self.v5 - self.x6) + v5dot - self.f2 )
        self.u = np.linalg.inv(self.jvec)* ( np.dot(np.transpose(self.gblock),self.varr1) + np.array([v4dot,v6dot]) + self.A7*self.varr2 )

    def integrate(self,q1,q2):
        return q2*self.dt + q1

    def derivative(self,q1,q2):
        return (q1-q2)/self.dt

    
    def compute(self):
        
        T = np.linspace(0,20,200)
        
        for i in range(1,200):
            self.t = T[i]
            self.v1o,self.v2o,self.v3o,self.v4o,self.v5o,self.v6o,self.uo = self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u
            self.update()
            self.calc_input()




            # self.t1 = T[i+1]
            # self.x0ddot, self.x0d = self.desiredVals(self.t0)
            # self.x1ddot, self.x1d = self.desiredVals(self.t1)
            # self.v1a = self.A1*(self.x0d - self.x1) + self.x0ddot
            # self.v1b = self.A1*(self.x1d - self.x1) + self.x1ddot
            # self.v1d = self.integrate(self.v1a, self.v1b)
            
            # self.v2a = np.linalg.inv(self.g0)*((self.x0d - self.x1) + self.A2*(self.v1 - self.x2) + self.v1d - self.f0)


            np.array([self.F1, self.F2, self.F3, self.F4]) = self.integrate(self.x7,self.x7dot)
            np.array([self.psidot, self.zdot]) = self.integrate(self.x6,self.x6dot)
            np.array([self.psi, self.z]) = self.integrate(self.x5,self.x5dot)
            np.array([self.thetadot, self.phidot]) = self.integrate(self.x4,self.x4dot)
            np.array([self.theta, self.phi]) = self.integrate(self.x3,self.x3dot)
            np.array([self.xdot, self.ydot]) = self.integrate(self.x2,self.x2dot)
            np.array([self.x,self.y]) = self.integrate(self.x1,self.x1dot)




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