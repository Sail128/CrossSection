import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

import numpy as np
import matplotlib.pyplot as plt

#using the easy boom analysis
class WinboxSection():
    def __init__(self):
        self.Nstringerstop = 30
        self.Nstringersbottom = 10
        self.Astringertop = 0.001 # [m^2]
        self.Astringerbottom = 0.001 # [m^2]
        self.Asparcap = 0.005 # [m^2]
        self.As_sparcaps = [0.005, #top right sparcap
                            0.003, #top left sparcap
                            0.003, #bottom left
                            0.005] #bottom right
        self.Hfront = 0.6
        self.Hback = 0.45
        self.tSpar = 0.02
        self.tskin = 0.001
        self.c = 3.0
        self.boomList = []
        self.tList = [] # list of the thickness of each section between boomlist[i] and boomlist[i-1]
        self.Alist = [] #  just a list of the stringer crossectional area with out panel addition
        self.centroid = []#if the list is empty than the centroid is not jet calculated
        self.Ixx = None #not calculated if None
        self.Iyy = None # not calculated is None
        self.Iyy = None # not calculated is None
        self.shearCenter = []
        self.NA = None # tuple(xs, ys, alpha)
        self.createGeometry()

    def createGeometry(self):
        # tho boom list is ordered clockwise starting at the top left corner
        # boom is defined as [x,y,A]
        self.boomList = [] # empty the boom list
        # open toppanel 
        self.boomList.append([0.0,self.Hfront/2,self.As_sparcaps[0]])     # top panel sparcap
        self.tList.append(self.tSpar)
        #calculate the top panel stringer
        Lseg = self.c/(self.Nstringerstop+1)
        def ytop(x):
            a = (self.Hback/2-self.Hfront/2)/(self.c-0.0)
            b = self.Hfront/2
            return a*x + b
        for i in range(self.Nstringerstop):
            x = (i+1)*Lseg
            self.boomList.append([x,ytop(x),self.Astringertop])
            self.tList.append(self.tskin)
        #close the top panel
        self.boomList.append([self.c,self.Hback/2,self.As_sparcaps[1]])   # top panel sparcap
        self.tList.append(self.tskin)
        # open bottom panel
        self.boomList.append([self.c,-self.Hback/2,self.As_sparcaps[2]])  # bottom panel sparcap
        self.tList.append(self.tSpar)
        #calculate the top panel stringer
        Lseg = self.c/(self.Nstringersbottom+1)
        def ybottom(x):
            a = (-self.Hback/2+self.Hfront/2)/(self.c-0.0)
            b = -self.Hfront/2
            return a*x + b
        for i in range(self.Nstringersbottom):
            x = self.c-(i+1)*Lseg
            self.boomList.append([x,ybottom(x),self.Astringerbottom])
            self.tList.append(self.tskin)
        # close bottom panel
        self.boomList.append([0.0,-self.Hfront/2,self.As_sparcaps[3]])    # bottom panel sparcap   
        self.tList.append(self.tskin)

        self.Alist = [b[2] for b in self.boomList]
        print(self.Alist)
        for i in range(20):
            a,b,c,NA =self.calculateParameters()
            xs,ys,A = NA
            print(A)
            Y = lambda x : A*x + self.shearCenter[1]

            for i in range(len(self.boomList)):
                b1 = self.boomList[i-1]
                b2 = self.boomList[i]
                b3 = self.boomList[(i+1)%len(self.boomList)]
                
                t12 = self.tList[i]
                t23 = self.tList[(i+1)%len(self.boomList)]
                S12 = np.sqrt((b2[1]-b1[1])**2+(b2[0]-b1[0])**2)
                S23 = np.sqrt((b3[1]-b2[1])**2+(b3[0]-b2[0])**2)
                self.boomList[i][2] = self.Alist[i]
                self.boomList[i][2] += t12*S12/6*(2+(b1[1]-Y(b1[0]))/(b2[1]-Y(b2[0])))
                self.boomList[i][2] += t23*S23/6*(2+(b1[1]-Y(b1[0]))/(b2[1]-Y(b2[0])))

    def getCentroid(self):
        if len(self.boomList) == 0: self.createGeometry()
        Atot = 0.0
        Xc = 0.0
        Yc = 0.0
        for boom in self.boomList:
            Atot +=boom[2]
            Xc += boom[0]*boom[2]
            Yc += boom[1]*boom[2]
        self.centroid = (Xc/Atot,Yc/Atot)
        return (Xc/Atot,Yc/Atot,Atot)

    def getI(self):
        if len(self.centroid) != 2:
            self.getCentroid()
        Ixx = 0.0
        Iyy = 0.0
        Ixy = 0.0
        for boom in self.boomList:
            I = np.pi/4 * (boom[2]/np.pi)**2
            dx = boom[0]-self.centroid[0]
            dy = boom[1]-self.centroid[1]
            #print(I, dy, boom[2]* dy**2)
            Ixx += I + boom[2]* dy**2
            Iyy += I + boom[2]* dx**2
            Ixy += boom[2]*dx*dy
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Ixy = Ixy
        return (Ixx, Iyy, Ixy)

    #returns a list if the shearflows between each boom, where flow i is between boom i and i-1
    def getUnitShearflow(self, Sx, Sy):
        if len(self.centroid) != 2:
            self.getCentroid()
        if self.Ixx ==None or self.Iyy == None or self.Ixy ==None:
            self.getI()
        a = (Sx*self.Ixx-Sy*self.Ixy)/(self.Ixx*self.Iyy-self.Ixy**2)
        b = (Sy*self.Iyy-Sx*self.Ixy)/(self.Ixx*self.Iyy-self.Ixy**2)
        xc = self.centroid[0]
        yc = self.centroid[1]
        Bx = 0.0
        By = 0.0
        for boom in self.boomList:
            Bx += boom[2]*(boom[0]-xc)
            By += boom[2]*(boom[1]-yc)
        qlist = []
        q_t = 0.0
        S_t = 0.0
        isbottom = 0
        for i in range(len(self.boomList)):
            b1 = self.boomList[i-1]
            b2 = self.boomList[i]
            x1 = b1[0]-xc
            y1 = b1[1]-xc
            S = np.sqrt((b2[1]-b1[1])**2+(b2[0]-b1[0])**2)
            try:
                alpha = (b2[1]-b1[1])/(b2[0]-b1[0])
                c = self.tList[i]*(0.5*np.cos(np.arctan(alpha))*S**2+x1*S)
                d = self.tList[i]*(0.5*np.sin(np.arctan(alpha))*S**2+y1*S)
                if isbottom == 2:
                    c = self.tList[i]*(0.5*-1*np.cos(np.arctan(alpha))*S**2+x1*S)
                    d = self.tList[i]*(0.5*-1*np.sin(np.arctan(alpha))*S**2+y1*S)
                c +=Bx
                d +=By
            except ZeroDivisionError:
                c = self.tList[i]*(x1*S)+Bx
                if isbottom == 0:
                    d = self.tList[i]*(0.5*S**2+y1*S)+By
                else:
                    d = self.tList[i]*(-0.5*S**2+y1*S)+By
                isbottom +=1

            qs = a*c-b*d
            qlist.append(qs)
            q_t += qs/self.tList[i]
            S_t += S/self.tList[i]

        return [qlist, q_t/S_t]

    def getShearCenter(self):
        qSy = self.getUnitShearflow(0,1)
        qSx = self.getUnitShearflow(1,0)
        #print(qSx)
        x_sc = 0.0
        y_sc = 0.0 #(qSx[0][0]+qSx[1])*(self.boomList[0][1]-self.boomList[-1][1])*(self.boomList[0][0]-self.centroid[0])
        for i in range(len(self.boomList)):
            b1 = self.boomList[i-1]
            b2 = self.boomList[i]
            S = np.sqrt((b2[1]-b1[1])**2+(b2[0]-b1[0])**2)
            x_sc += (qSy[0][i]+qSy[1])*(S)*((b1[0]+b2[0])/2-self.centroid[0])
            y_sc += (qSx[0][i]+qSx[1])*(S)*((b1[1]+b2[1])/2-self.centroid[1])
        #print(x_sc, y_sc)
        self.shearCenter = (self.centroid[0] + x_sc, self.centroid[1]+y_sc)
        return self.shearCenter


    def getNA(self):
        if len(self.shearCenter) == 0:
            self.getShearCenter()

        A = -(self.Ixy/self.Iyy)
        xs = np.arange(-3, self.c+3)
        ys = xs*A + self.shearCenter[1]
        self.NA = (xs, ys, A)
        return self.NA

    def plotGeometry(self):
        xs = [i[0] for i in self.boomList]
        ys = [i[1] for i in self.boomList]
        As = [i[2]*10000 for i in self.boomList]
        plotx = xs
        ploty = ys
        plotx.append(xs[0])
        ploty.append(ys[0])
        plt.plot(plotx,ploty, color="red")
        plt.scatter(xs,ys,s=As)

        if len(self.centroid) == 2:
            plt.scatter(self.centroid[0],self.centroid[1], marker="+" ,s = 100,label="centroid")
        if len(self.shearCenter) == 2:
            plt.scatter(self.shearCenter[0],self.shearCenter[1], marker="x" ,s = 100,label="shear center")
        if self.NA != None:
            xs,ys, alpha = self.NA
            plt.plot(xs,ys, linestyle="--", label="Neutral axis at {:0.2e}[deg]".format(np.arctan(alpha)))
        plt.axis("equal")
        plt.xlim( self.shearCenter[0]-0.1 if self.shearCenter[0]<0 else -0.1, self.c+0.1)
        plt.legend()
        plt.show()
        plt.close()

    def calculateParameters(self):
        centroid = self.getCentroid()
        inertiaMoments = self.getI()
        SC = self.getShearCenter()
        NA = self.getNA()
        return (centroid, inertiaMoments, SC ,NA)

if __name__ == "__main__":
    test = WinboxSection()
    test.createGeometry()
    c,I,SC,NA = test.calculateParameters()
    Ixx,Iyy,Ixy = I
    c_x,c_y,A = c
    print("centroid: ",(c_x,c_y))
    print("crossection: ",A)
    print("Ixx: ", Ixx)
    print("Iyy: ", Iyy)
    print("Ixy: ", Ixy)
    print("shear center: ", SC)
    print("d/dx(NA): ", NA[2])
    test.plotGeometry()