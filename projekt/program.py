from transformacje import Transformacje

print('''Wybierz operację na współrzędnych:\n
      1 - ze wsp geocentrycznych do wsp geodezyjnych
      2 - ze wsp geodezyjnych do wsp geocentrycznych
      3 - do wsp topocentrycznych (NEU)
      4 - do układu 2000
      5 - do układu 1992
      6 - odległosc 2D
      7 - odległosc 3D ''')
      
opcje = input()

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    
    if '1' in opcje:
        print("Wprowadź współrzędne do transformacji (X,Y,Z)")
        X = input('X: ')
        Y = input('Y: ')
        Z = input('Z: ')
        phi, lam, h = geo.xyz2flh(X, Y, Z)
        print(phi, lam, h)
        
    if '2' in opcje:
        print("Wprowadź współrzędne do transformacji (phi,lam,h)")
        phi = input('phi: ')
        lam = input('lambda: ')
        h = input('h: ')
        X, Y, Z = geo.flh2xyz(phi,lam,h)
        print(X, Y, Z)
        
    if '3' in opcje:
        print("Wprowadź współrzędne do transformacji (X,Y,Z)")
        X = input('X: ')
        Y = input('Y: ')
        Z = input('Z: ')
        n,e,u = geo.xyz2neu(X,Y,Z)
        print(n,e,u)
        
    if '4' in opcje:
        print("Wprowadź współrzędne do transformacji (phi,lam)")
        phi = input('phi: ')
        lam = input('lambda: ')
        x2000,y2000 = geo.u2000(phi,lam)
        print(x2000,y2000)
        
    if '5' in opcje:
        print("Wprowadź współrzędne do transformacji (phi,lam)")
        phi = input('phi: ')
        lam = input('lambda: ')
        x1992,y1992 = geo.u1992(phi,lam)
        print(x1992,y1992)
        
    if '6' in opcje:
        print("Wprowadź współrzędne do transformacji (X,Y,X2,Y2)")
        X = input('X: ')
        Y = input('Y: ')
        X2 = input('X2: ')
        Y2 = input('Y2: ')
        odl2D = geo.odl2D(X, Y, X2, Y2)
        print(odl2D)
        
    if '7' in opcje:
        print("Wprowadź współrzędne do transformacji (X,Y,Z,X2,Y2,Z2)")
        X = input('X: ')
        Y = input('Y: ')
        Z = input('Z: ')
        X2 = input('X2: ')
        Y2 = input('Y2: ')
        Z2 = input('Z2: ')
        odl3D = geo.odl3D(X,Y,Z,X2,Y2,Z2)
        print(odl3D)
