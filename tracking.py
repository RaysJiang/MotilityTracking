"""
0   Timepoint
1   Time [s]
2   Object ID
3   X
4   Y
5   Tracked Nuclei - Age [s]
6   Tracked Nuclei - Current Displacement X [um]
7   Tracked Nuclei - Current Displacement Y [um]
8   Tracked Nuclei - Track Point X
9   Tracked Nuclei - Track Point Y
10  Tracked Nuclei - Current Step Size [um]
11  Tracked Nuclei - Current Speed [um/s]
...
-2  Tracks: Tracked Nuclei - Accumulated Distance [um]
-1  Tracks: Tracked Nuclei - Displacement [um]
"""


import math
import copy

class track():
    def __init__(self):
            self.uXY = {}
            self.uS = {}
            self.uTracking = {}
            self.uAngle = {}
            self.uDistance = {}
            self.IDs = []

    def LoadInput(self,InputFile):
        F = open(InputFile,'r')
        f = F.read()
        f = f.replace('\r','\n')
        f = f.split('\n')
    
        n = 0
        uXY = self.uXY
        uS = self.uS
        uDistance = self.uDistance
        for i in f:
            i = i.strip()
            if i != '':
                n += 1
            if n > 1:
                x = i.split('\t')
                T = int(x[0])
                ID = int(x[2])
                X = int(x[3])
                Y = int(x[4])
                S = float(x[11])    
                Dis = float(x[-1])
            
                if uXY.has_key(ID):
                    tmp = uXY[ID]
                    tmp[T]=(X,Y)
                else:
                    uXY[ID]={}
                    tmp = uXY[ID]
                    tmp[T]=(X,Y)
                
                if uS.has_key(ID):
                    tmp = uS[ID]
                    tmp.append(S)
                else:
                    uS[ID]=[]
                    tmp = uS[ID]
                    tmp.append(S)  
                    
                if uDistance.has_key(ID):
                    pass
                else:
                    uDistance[ID] = Dis  
                             
                
        F.close()
    
        IDs = uXY.keys()
        IDs.sort()
        self.IDs = IDs
        for i in IDs:
            utmp = uXY[i]
            utmp_centered = self.recenter(utmp)
            uXY[i] = utmp_centered
            print i,utmp,'\n',i,utmp_centered,'\n\n\n'
    
        print n, 'lines parsed'
        print len(self.IDs),'objects'

        
    def tracking(self):
        uXY = self.uXY
        uS = self.uS
        uDistance = self.uDistance
        uAngle = self.uAngle
        uTracking = self.uTracking
        IDs = self.IDs

        for i in IDs:
            S = uS[i]
        
            uXYtmp = uXY[i]
            uPoints = self.ThreePoints(uXYtmp)
            vAngles = self.angles(uPoints)
        
            uAngle[i]=vAngles
        
            Smean = average(S)
            Amean = average(vAngles)
            Dis = uDistance[i]
            T = len(S)
        
            uTracking[i]=[Smean,Amean,Dis,T]
            #print i,Smean,Amean
            
            
        G = open('Output_tracking_summary.txt','w')
        G.write('SpeedRange(um),Count,Percentage\n'.replace(',','\t'))
        Sall = [n[0] for n in uTracking.values()]   
        H = histogram(Sall)  
        Bins = H.keys()
        Bins.sort()
        N = float(len(IDs))
        for bin in Bins:
            a,b = bin
            a,b= str(round(a,3)),str(round(b,3))
            C = H[bin]
            G.write('%s_%s\t%d\t%f\n'%(a,b,C,C*100/N))
            print a,b,C
            
       
        G.write('\n\n\n'+'#'*30+'\n\n\n')  
        G.write('TrackingID,TimePoints,SpeedMean,AngleMean,TotalDistance\n'.replace(',','\t'))
        for i in IDs:
            Smean,Amean,Dis,T = uTracking[i]
            G.write('%d\t%d\t%f\t%f\t%f\n'%(i,T,Smean,Amean,Dis))
        
        
        G.close()
        
 
    
    
    def selecting(self,speedlimit):
        uXY = self.uXY
        uS = self.uS
        uAngle = self.uAngle
        uTracking = self.uTracking
        IDs = self.IDs
        

        G = open('output_SelectedSporozoites.txt','w')
        G.write('ID,Time,X,Y,Speed,Angle\n'.replace(',','\t'))
        c = 0
        for i in IDs:
            S = uS[i]
            A = ['NA']*3+uAngle[i]
            Smean,Amean,Dis,T = uTracking[i]
            if Smean >= speedlimit:
                c += 1
                uXYtmp = uXY[i]
                T = uXYtmp.keys()
                T.sort()
                n = 0
                for t in T:
                    x,y = uXYtmp[t]
                    s = S[n]
                    a = A[n]
                    n += 1
                    tmpoutput = strlist((i,t,x,y,s,a))
                    G.write('\t'.join(tmpoutput)+'\n')
        G.close()  
        print c,'selected'
        print len(IDs),'in total'          

     
    
    def recenter(self,u):
        tmp = u[0]
        X0,Y0 = tmp
        X0 = 0-X0
        Y0 = 0-Y0
    
        uCenter = {}
        Time = u.keys()
        Time.sort()
        for t in Time:
            X,Y = u[t]
            uCenter[t]=[X+X0,Y+Y0]
        
        return uCenter  
    
    def ThreePoints(self,u):
        Time = u.keys()
        Time.sort()
        n = 0
        u3points = {}
        points = []
        for t in Time:
            n += 1
            points.append(u[t])
            if n >= 3:
                u3points[t]=copy.copy(points)
                del points[0]
        return u3points 
    
    def triangle(self,P1,P2,P3):
        P1x,P1y = P1
        P2x,P2y = P2
        P3x,P3y = P3
        #sqrt((P1x - P2x)2 + (P1y - P2y)2)
        D12 = math.sqrt(math.pow(P1x-P2x,2)+math.pow(P1y-P2y,2))
        D23 = math.sqrt(math.pow(P2x-P3x,2)+math.pow(P2y-P3y,2))
        D13 = math.sqrt(math.pow(P1x-P3x,2)+math.pow(P1y-P3y,2))
        #=(P122 + P132 - P232) / (2 * P12 * P13)
        tri_beta = (math.pow(D12,2)+math.pow(D13,2)-math.pow(D23,2))/(2*D12*D13)
        tri_alpha = (math.pow(D23,2)+math.pow(D13,2)-math.pow(D12,2))/(2*D23*D13)
        # arcos
        beta = math.acos(tri_beta)
        alpha = math.acos(tri_alpha)
        #degrees
        beta = math.degrees(beta)
        alpha = math.degrees(alpha)
        gamma = 180-beta-alpha
    
        return gamma               
        
    def angles(self,u):
        vAngles = []
        Time = u.keys()
        Time.sort()
        for t in Time:
            points = u[t]
            P1,P2,P3 = points
            gamma = 0
            try:
                gamma = self.triangle(P1,P2,P3)    
            except :
                pass
                #print points,'error!'   
            
            vAngles.append(gamma)
        return vAngles       
        
        
         
def strlist(x):
    y = []
    for i in x:
        y.append(str(i))
    return y        
    
def average(X):
    S = sum(X)
    L = float(len(X))
    mean = None
    if L != 0:
        mean = S/L
    else:
        pass    
    return mean
    
def histogram(x):
    bins = map(lambda x: x/10.0, range(0, 35, 1))
    A = bins[:-1]
    B = bins[1:]
    u = {}
    n = 0
    for a in A:
        b = B[n]
        u[(a,b)]= 0
        n += 1
    
    bins = u.keys()
    bins.sort()
    bin_last = bins[-1]
    bin_limit = bin_last[1]    
    for i in x:
        for bin in bins:
            a,b = bin
            if i >= a and i < b:
                u[bin]= u[bin]+1
            elif i >= bin_limit:
                u[bin_last] = u[bin_last] + 1    
    
    return u            
            

        
    
if __name__ == '__main__':
    InputFile= 'Input_PointsTracking.txt'
    speedlimit = 0.5
    T = track()
    T.LoadInput(InputFile)   
    T.tracking() 
    T.selecting(speedlimit)
                    
                        