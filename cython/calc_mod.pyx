def calc(img,frecuencia,espacio):
    h,c,w=img.shape #h para filas c para columnas w para numero de matrices 1 con cada color
    ima=[0]*h
    for i in range(len(ima)):
        ima[i]=[0]*c
    for i in range(h):
        for j in range(c):
            promedio=(int(img[i][j][0])+int(img[i][j][1])+int(img[i][j][2]))/3
            ima[i][j]=int(promedio) #valor de gris igual al promedio de RGB
    imb=[None]*h
    for i in range(len(imb)):
        imb[i]=[None]*c
    if ima[0][0]<128: #binarizando para fondo negro
        for i in range(h):
            for j in range(c):
                if ima[i][j]<80:
                    imb[i][j]=0
                else:
                    imb[i][j]=255
    else: #binarizando para fondo blanco
        for i in range(h):
            for j in range(c):
                if ima[i][j]>200:
                    imb[i][j]=0
                else:
                    imb[i][j]=255
    def label(img):#colorea cada objeto
        i=1
        intervalos=[]
        while i<h:
            if 255 not in img[i] and 255 in img[i-1]:
                intervalos.append(i)
            i+=1
        n=1
        intervalos=intervalos
        while n<len(intervalos):
            k=intervalos[n-1]
            while k<intervalos[n]:
                for j in range(c):
                    if img[k][j]!=0:
                        img[k][j]=n
                k+=1
            n+=1
        return img,n-1,intervalos
    def intervalos_reversa(img):
        i=h-2
        intervalos=[]
        while i>0:
            if 255 not in img[i] and 255 in img[i+1]:
                intervalos.append(i)
            i-=1
        intervalos.reverse()
        return intervalos
    def ccm(intervalos,ir):
        ccm=[]
        del intervalos[0]
        i=0
        while i<len(intervalos):
            ccm.append((int(intervalos[i]+ir[i]))/2)
            i+=1
        del ccm[0]
        return ccm
    ir=intervalos_reversa(imb)
    imc,n,i=label(imb) #el elemento 0 es la matriz, el 1 es el numero de objetos contados, el 2 el intervalo
    lc=ccm(i,ir)#lista coordenadas
    lc.remove(min(lc)) #quitando la primera bolita
    def mult(a,b):
        try:
            C = [[0 for fila in range(len(b[0]))] for columna in range(len(a))]
            for i in range(len(a)):
                for j in range(len(b[0])):
                    for k in range(len(a[0])):
                        C[i][j]+=a[i][k]*b[k][j]
            return C
        except: #por si la matriz b es un vector nx1, ya que len(float) es error
            C = [0]*len(a)
            for i in range(len(a)):
                for k in range(len(a[0])):
                    C[i]+=a[i][k]*b[k]
            return C
    def transpuesta(m):
        return [list(i) for i in zip(*m)]
    def menor(m,i,j): #cuando uno tapa la fila i esima con la columna j esima
        return [fila[:j] + fila[j+1:] for fila in (m[:i]+m[i+1:])]
    def determinante(m):
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]
        determinant = 0
        for c in range(len(m)):
            determinant += ((-1)**c)*m[0][c]*determinante(menor(m,0,c))
        return determinant
    def inversa(m):
        determinant = determinante(m)
        if len(m) == 2:
            return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]
        #matriz de cofactores
        cofactores = []
        for r in range(len(m)):
            filacofactores = []
            for c in range(len(m)):
                menorr = menor(m,r,c)
                filacofactores.append(((-1)**(r+c)) * determinante(menorr))
            cofactores.append(filacofactores)
        cofactores = transpuesta(cofactores)
        for r in range(len(cofactores)):
            for c in range(len(cofactores)):
                cofactores[r][c] = cofactores[r][c]/determinant
        return cofactores
    def minimos_cuadrados(t,y):
        Xt=[[],[],[]]
        for i in range(len(t)):
            Xt[0].append(1)
            Xt[1].append(t[i])
            Xt[2].append(t[i]**2)
        X=[list(i) for i in zip(*Xt)]
        a = mult(mult(inversa(mult(Xt,X)),Xt),y)    
        return a
    tiempos=list(range(1,len(lc)+1))
    coeficientes=minimos_cuadrados(tiempos,lc) #y(t)=c+bt+0.5at2 lc son la lista de las posiciones
    gravedad=coeficientes[-1]*2.0
    def SI(g,Uespacio,freq): #Uespacio es el tamano de cada pixel
        return g*(freq**2)*Uespacio #para que de en mm/s^2 si las entradas estan en mm y s
    return SI(gravedad,espacio,frecuencia)
    
