import matplotlib.pyplot as plt
import numpy as np
from numpy import sin as S 
from numpy import cos as C 
from numpy import tan as T
from numpy.core.shape_base import block
from numpy.lib.index_tricks import ix_ 
from scipy.linalg import block_diag
import time

class controller:
    def __init__(self):
        self.g = 9.81
        self.m = 2.0
        self.It = np.diag([1.2416, 1.2416,2.4832])
        self.Ix, self.Iy, self.Iz = 1.2416, 1.2416, 2.4832
        self.F1, self.F2, self.F3, self.F4 = self.m*self.g/4,self.m*self.g/4,self.m*self.g/4,self.m*self.g/4
        self.F1dot, self.F2dot, self.F3dot, self.F4dot = 0.0,0.0,0.0,0.0
        
        self.d = 0.2
        self.c = 0.01
        
        self.kt = np.diag([0.01, 0.01, 0.01])
        self.kr = np.diag([0.001, 0.001, 0.001])
        self.A1,self.A2,self.A3,self.A4,self.A5,self.A6 = np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),np.diag([3,3]),
        self.A7 = np.diag([3,3,3,3])
        self.dt = 0.1
        self.t = 0.0

        self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u = np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0],[0.0],[0.0]])
        # self.fx, self.fy, self.fz = 0.0, 0.0, -self.g
        self.farr1 = np.array([[0],[0],[-self.g]])

        self.neta = np.array([[0],[0],[0]])
        self.netadot = np.array([[0],[0],[0]])
        self.eta = np.array([[0],[0],[0]])
        self.etadot = np.array([[0],[0],[0]])

        self.x1 = np.array([[0],[0]])
        # print(self.x1)
        self.x2 = np.array([[0], [0]])
        self.x3 = np.array([[0], [0]])
        self.x4 = np.array([[0], [0]])
        self.x5 = np.array([[0], [0]])
        self.x6 = np.array([[0], [0]])
        self.x7 = np.array([[0], [0], [0], [0]])
        # self.u = np.array([[0], [0], [0], [0]])


    def paramters_update(self):
        self.g0 = ((self.F1+self.F2+self.F3+self.F4)/self.m)*np.array([[S(self.x5[0]), C(self.x5[0])],[-C(self.x5[0]), S(self.x5[0])]])
        # print(self.g0)
        self.g1 = np.array([[1/self.Ix , S(self.x3[0])*T(self.x3[1])/self.Iy],[0, C(self.x3[0])/self.Iy]])
        self.g2 = np.array([[C(self.x3[0])/(self.Iz*C(self.x3[1])) , 0],[0, C(self.x3[0])*C(self.x3[1])/self.m]])

        self.pivec0 = np.array([[S(self.x3[0])],[C(self.x3[0])*S(self.x3[1])]])
        self.pivec1 = np.array([[self.d*(self.F2-self.F4)],[self.d*(self.F3-self.F1)]])
        self.pivec2 = np.array([[self.c*(self.F1-self.F2+self.F3-self.F4)],[self.F1+self.F2+self.F3+self.F4]])

        self.Rt = np.array([[C(self.x3[0])*C(self.x5[0]), (S(self.x3[0])*S(self.x3[1])*C(self.x5[0]) - C(self.x3[0])*S(self.x5[0])), (C(self.x3[0])*S(self.x3[1])*C(self.x5[0]) + S(self.x3[0])*S(self.x5[0]))], [C(self.x3[1])*S(self.x5[0]), (S(self.x3[0])*S(self.x3[1])*S(self.x5[0]) + C(self.x3[0])*C(self.x5[0])), (C(self.x3[0])*S(self.x3[1])*S(self.x5[0]) - S(self.x3[0])*C(self.x5[0]))], [-S(self.x3[0]), S(self.x3[0])*C(self.x3[1]), C(self.x3[0])*C(self.x3[1])]])
        self.Rr = np.array([[1,0,-S(self.x3[1])], [0,C(self.x3[0]),C(self.x3[1])*S(self.x3[0])], [0,-S(self.x3[0]),C(self.x3[0])*C(self.x3[1])]])
        self.dRrbydphi = np.array([[0,0,0], [0,-S(self.x3[0]),C(self.x3[1])*C(self.x3[0])], [0,-C(self.x3[0]),-S(self.x3[0])*C(self.x3[1])]])
        self.dRrbydtheta = np.array([[0,0,-C(self.x3[1])], [0,0,-S(self.x3[1])*S(self.x3[0])], [0,0,-C(self.x3[0])*S(self.x3[1])]])

        # print(np.shape(self.etadot), "etadot")
        # print(np.shape(self.Rr), "Rr")
        # print(np.shape(self.It), "It")
        # print(np.shape(self.Rr.dot(self.etadot)))
        # print(np.shape(self.It.dot(self.Rr).dot(self.etadot)))
        
        self.farr2 = -np.linalg.inv(np.dot(self.It,self.Rr)).dot((self.It.dot(self.x4[0]*self.dRrbydphi+self.x4[1]*self.dRrbydtheta).dot(self.etadot) - np.cross(self.Rr.dot(self.etadot), self.It.dot(self.Rr).dot(self.etadot),axis=0))) + np.array([[(self.c/self.Iz)*C(self.x3[0])*T(self.x3[1])*(self.F1-self.F2+self.F3-self.F4)],[(-self.c/self.Iz)*S(self.x3[0])*(self.F1-self.F2+self.F3-self.F4)],[(self.d/self.Iy)*(S(self.x3[0])/C(self.x3[1]))*(self.F3-self.F1)]])
        # print(self.farr2)
        # print((self.It.dot(self.x4[0]*self.dRrbydphi+self.x4[1]*self.dRrbydtheta).dot(self.etadot)))# - 
        # print(np.cross(self.Rr.dot(self.etadot), self.It.dot(self.Rr).dot(self.etadot)))
        # print(np.array([[(self.c/self.Iz)*C(self.x3[0])*T(self.x3[1])*(self.F1-self.F2+self.F3-self.F4)],[(-self.c/self.Iz)*S(self.x3[0])*(self.F1-self.F2+self.F3-self.F4)],[(self.d/self.Iy)*(S(self.x3[0])/C(self.x3[1]))*(self.F3-self.F1)]]))
        # self.farr2 = np.array([[self.farr2[0]],[self.farr2[1]],[self.farr2[2]]])
        # print(self.farr2[1])
        # print(self.farr2[0])
        # print(self.farr2)


        self.f0 = np.vstack((self.fx,self.fy))
        self.f1 = np.vstack((self.farr2[0],self.farr2[1]))
        self.f2 = np.vstack((self.farr2[2],self.fz))
        # print(self.f1)
        self.J0 = np.array([[C(self.x3[0]),0],[-S(self.x3[0])*S(self.x3[1]), C(self.x3[0])*C(self.x3[1])]])
        self.J1 = np.array([[0,self.d,0,-self.d],[-self.d,0,self.d,0]])
        self.J2 = np.array([[self.c,-self.c,self.c,-self.c],[1,1,1,1]])

        self.jvec = np.vstack((self.J1,self.J2))
        # print(np.shape(self.jvec))
        # print(np.shape(self.J1), np.shape(self.J2))

        self.gblock = block_diag(self.g1,self.g2)

        self.x1dot = self.x2
        self.x2dot = self.f0 + self.g0.dot(self.pivec0)
        self.x3dot = self.x4
        self.x4dot = self.f1 + self.g1.dot(self.pivec1)
        self.x5dot = self.x6
        self.x6dot = self.f2 + self.g2.dot(self.pivec2)
        self.x7dot = self.u

        # print(self.gblock)

    def integrate(self,q1,q2):
        return q2*self.dt + q1

    def derivative(self,q1,q2):
        # print((q1-q2)/self.dt)
        return (q1-q2)/self.dt

    def calculate(self):
        self.v1 = self.A1@((self.x1d-self.x1)) + self.x1ddot
        v1dot = self.derivative(self.v1,self.v1o)
        self.v2 = np.linalg.inv(self.g0)@(((self.x1d - self.x1)) + self.A2@((self.v1 - self.x2)) + v1dot - self.f0 )
        v2dot = self.derivative(self.v2,self.v2o)
        self.v3 = np.linalg.inv(self.J0)@(np.transpose(self.g0)@((self.v1 - self.x2)) + self.A3@(self.v2 - self.pivec0) + v2dot )
        v3dot = self.derivative(self.v3,self.v3o)
        self.v4 = np.linalg.inv(self.g1)@( np.transpose(self.J0)@(self.v2 - self.pivec0) + self.A4@(self.v3 - self.x4) + v3dot - self.f1 )
        v4dot = self.derivative(self.v4,self.v4o)
        self.v5 = self.A5@(self.x5d - self.x5) + self.x5ddot
        v5dot = self.derivative(self.v5,self.v5o)
        self.varr1 = np.vstack((self.v3 - self.x4, self.v5 - self.x6))
        self.v6 = np.linalg.inv(self.g2)@( (self.x5d - self.x5) + self.A6@(self.v5 - self.x6) + v5dot - self.f2 )
        v6dot = self.derivative(self.v6,self.v6o)
        self.varr2 = np.vstack((self.v4 - self.pivec1,self.v6 - self.pivec2))
        self.u = np.linalg.inv(self.jvec).dot(((self.varr1.dot(np.transpose(self.gblock))) + np.vstack((v4dot,v6dot)) + self.A7.dot(self.varr2)))
        self.x7dot = self.u

    def mainloop(self):
        T = np.linspace(0,20,200)        
        plt.figure(figsize=(10,5))
        ax = plt.axes(projection ='3d')

        for i in range(1,200):
            self.t = T[i]
            self.v1o,self.v2o,self.v3o,self.v4o,self.v5o,self.v6o,self.uo = self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u

            # DESIRED TRAJECTORY PARAMETERS
            self.x5d = np.array([[0], [10*S(self.t)]])
            self.x5ddot = np.array([[0], [10*C(self.t)]])
            self.x1ddot, self.x1d = np.array([[0], [0]]), np.array([[0], [0]])
            self.paramters_update()

            self.F1, self.F2, self.F3, self.F4 = self.integrate(self.x7,self.x7dot)[:,0]
            print(self.integrate(self.x7,self.x7dot)[:,0])
            self.x6[0], self.x6[1] = self.integrate(self.x6,self.x6dot)[:,0]
            self.x5[0], self.x5[1] = self.integrate(self.x5,self.x5dot)[:,0]
            self.x4[0], self.x4[1] = self.integrate(self.x4,self.x4dot)[:,0]
            self.x3[0], self.x3[1] = self.integrate(self.x3,self.x3dot)[:,0]
            self.x2[0], self.x2[1] = self.integrate(self.x2,self.x2dot)[:,0]
            self.x1[0],self.x1[1] = self.integrate(self.x1,self.x1dot)[:,0]
            self.calculate()

            ax.plot(self.x, self.y, self.z, c='lightblue',marker='o')
            ax.plot(0, 0, 10*S(self.t), c='red',marker='o')
            plt.show(block=False)
            plt.pause(0.01)
            # time.sleep(0.1)
        plt.show()

if __name__ == '__main__':
    nlc = controller()
    nlc.mainloop()