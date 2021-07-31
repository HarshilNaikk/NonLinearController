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
