from math import sin, cos, sqrt, atan, degrees, tan, radians, pi
import numpy as np


class Transformacje:
    def __init__(self, model: str = "wgs84"):
        if model == "wgs84":
            self.a = 6378137
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
        
    def Np(self, phi):
        N = self.a/(1 - self.ecc2*(sin(phi))**2)**(0.5)
        return N

    def xyz2flh(self, X, Y, Z):
        """
        Funkcja przekształcająca współrzędne geocentryczne (X,Y,Z) we współrzędne geodezyjne (phi, lambda, h).
        
        Parameters
        ----------
        X
        Y
        Z

        Returns
        -------
        phi
        lambda
        h
        """
        X = float(X)
        Y = float(Y)
        Z = float(Z) 
        r = sqrt(X * 2 + Y * 2)
        fi_n = atan(Z / (r * (1 - self.ecc2)))
        eps = 0.000001 / 3600 * pi / 180
        fi = fi_n * 2
        while abs(fi_n - fi) > eps:
            fi = fi_n
            N = self.Np(fi_n)
            h = r / np.cos(fi_n) - N
            fi_n = atan(Z / (r * (1 - self.ecc2 * (N / (N + h)))))
        lam = atan(Y / X)
        h = r / cos(fi_n) - N
        return (f"{degrees(fi):.6f}", f"{degrees(lam):.6f}", f"{h:.3f}")
 
        
    def flh2xyz(self,phi, lam, hel):
        """
        Funkcja przekształcająca współrzędne geodezyjne (phi, lambda, h) we współrzędne geocentryczne (X,Y,Z).
        
        Parameters
        ----------
        phi
        lambda
        h elipsoidalne

        Returns
        -------
        X
        Y
        Z
        """
        phi = float(phi)
        lam = float(lam)
        hel = float(hel) 
        N = self.a/(1 - self.ecc2*(sin(phi))**2)**(0.5)
        X = (N+hel)*cos(phi)*cos(lam)
        Y = (N+hel)*cos(phi)*sin(lam)
        Z = (N*(1-self.ecc2)+hel)*sin(phi)
        return(f"{X:.3f}", f"{Y:.3f}", f"{Z:.3f}")
    
    def xyz2neu(self, X,Y,Z):
        """
        Funkcja wyznaczająca współrzędne topocentryczne (N,E,U).

        Parameters
        ----------
        X
        Y
        Z
       
        Returns
        -------
        N
        E
        U
        """
        X = float(X)
        Y = float(Y)
        Z = float(Z) 
        phi,lam,h = self.xyz2flh(X, Y, Z)
        phi = radians(phi)
        lam = radians(lam)
        R = np.array([(-sin(phi)*cos(lam), -sin(lam), cos(phi) * cos(lam)), (-sin(phi)*sin(lam), cos(lam), cos(phi)*sin(lam)), (cos(phi), 0, sin(phi))])
        n = R[:,0]
        e = R[:,1]
        u = R[:,2]
        return(n, e, u)
    
    def gauss_kruger(self, phi, lam):
        """
        Funkcja wyznaczająca współrzędne Gaussa-Krugera.
    
        """
        phi = float(phi)
        lam = float(lam)
        b2 = (self.a**2)*(1-self.ecc2)
        ep2 = ((self.a**2)-(b2))/(b2)
        t = tan(phi)
        n2 = ep2 * ((cos(phi))**2)
        N1 = self.a/(1 - self.ecc2*(sin(phi))**2)**(0.5)
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)
        A2 = (3/8) * (self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128))
        A4 = (15/256) * ((self.ecc2**2)+((3*(self.ecc2**3))/4))
        A6 = (35*(self.ecc2**3))/3072
        sig = self.a*(A0*phi - A2*sin(2*phi) + A4*sin(4*phi) - A6*sin(6*phi))
        if lam < 0.261799:
            L0 = 15
            nrS = 5
        elif lam > 0.261799 and lam < 0.314159:
            L0 = 18
            nrS = 6
        elif lam > 0.314159 and lam < 0.366519:
            L0 = 21
            nrS = 7
        elif lam > 0.366519:
            L0 = 24
            nrS = 8
        dL = lam - L0;
        xgk = sig + ((dL**2)/2)*N1*sin(phi)*cos(phi)*(1+((dL**2)/12)*((cos(phi))**2)*(5-t**2+9*n2+4*(n2**2))+((dL**4)/360)*((cos(phi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk = dL * N1 * cos(phi) * (1+((dL**2)/6)*((cos(phi))**2)*(1-t**2+n2)+((dL**4)/120)*((cos(phi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
        return(xgk, ygk, nrS)
    
    def u2000(self, phi,lam):
        """
        Funckja wyznaczająca współrzędne w układzie 2000.
        """
        phi = float(phi)
        lam = float(lam)
        xgk,ygk,nrS = self.gauss_kruger(phi, lam)
        x2000 = xgk * 0.999923
        y2000 = ygk * 0.999923 + ((nrS * 1000000)+500000)
        
        return(x2000, y2000)
    
    def u1992(self, phi, lam):
        """
        Funckja wyznaczająca współrzędne w układzie 1992.
        """
        phi = float(phi)
        lam = float(lam)
        xgk,ygk,nrS = self.gauss_kruger(phi, lam)
        x92 = (xgk * 0.9993) - 5300000
        y92 = (ygk * 0.9993) + 500000
        
        return(x92, y92)
    
    def odl3D(self,X,Y,Z,X2,Y2,Z2):
        """
        Obliczenie odległosci 3D.
        """
        X = float(X)
        Y = float(Y)
        Z = float(Z) 
        X2 = float(X2)
        Y2 = float(Y2)
        Z2 = float(Z2) 
        A = [X, Y, Z]
        B = [X2, Y2, Z2]
        od = sqrt( (A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2 )
        return(od)
    
    def odl2D(self,X,Y,X2,Y2):
        """
        Obliczenie odległosci 2D.
        """
        X = float(X)
        Y = float(Y)
        X2 = float(X2)
        Y2 = float(Y2)
        A = [X, Y]
        B = [X2,Y2]
        odl2d = sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2)
        return(odl2d)
    


    
        
    
