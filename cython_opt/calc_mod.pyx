def calc(img,frecuencia,espacio):
    cdef int h=img.shape[0]
    cdef int c=img.shape[1]
    cdef int w=img.shape[2]#h para filas c para columnas w para numero de matrices 1 con cada color
    cdef list ima=[0]*h
    cdef int i
    for i in range(len(ima)):
        ima[i]=[0]*c
    cdef int i
    for i in range(h):
        cdef int j
        for j in range(c):
            cdef int promedio=(int(img[i][j][0])+int(img[i][j][1])+int(img[i][j][2]))/3
            ima[i][j]=int(promedio) #valor de gris igual al promedio de RGB
    cdef list imb=[None]*h
    cdef int i
    for i in range(len(imb)):
        imb[i]=[None]*c
    if ima[0][0]<128: #binarizando para fondo negro
        cdef int i
        for i in range(h):
            cdef int j
            for j in range(c):
                if ima[i][j]<80:
                    imb[i][j]=0
                else:
                    imb[i][j]=255
    else: #binarizando para fondo blanco
        cdef int i
        for i in range(h):
            cdef int j
            for j in range(c):
                if ima[i][j]>200:
                    imb[i][j]=0
                else:
                    imb[i][j]=255
    def label(img):#colorea cada objeto
        cdef int i=1
        cdef list intervalos=[]
        while i<h:
            if 255 not in img[i] and 255 in img[i-1]:
                intervalos.append(i)
            i+=1
        cdef int n=1
        while n<len(intervalos):
            cdef int k=intervalos[n-1]
            while k<intervalos[n]:
                cdef int j
                for j in range(c):
                    if img[k][j]!=0:
                        img[k][j]=n
                k+=1
            n+=1
        return img,n-1,intervalos
    def intervalos_reversa(img):
        cdef int i=h-2
        cdef list intervalos=[]
        while i>0:
            if 255 not in img[i] and 255 in img[i+1]:
                intervalos.append(i)
            i-=1
        intervalos.reverse()
        return intervalos
    def ccm(list intervalos,list ir):
        cdef list ccm=[]
        del intervalos[0]
        cdef int i=0
        while i<len(intervalos):
            ccm.append((int(intervalos[i]+ir[i]))/2)
            i+=1
        del ccm[0]
        return ccm
    cdef list ir=intervalos_reversa(imb)
    imc,n,i=label(imb) #el elemento 0 es la matriz, el 1 es el numero de objetos contados, el 2 el intervalo
    cdef list lc=ccm(i,ir)#lista coordenadas
    lc.remove(min(lc)) #quitando la primera bolita
    cpdef list mult(list a,list b):
        try:
            cdef list C = [[0 for fila in range(len(b[0]))] for columna in range(len(a))]
            cdef int i
            for i in range(len(a)):
                cdef int j
                for j in range(len(b[0])):
                    cdef int k
                    for k in range(len(a[0])):
                        C[i][j]+=a[i][k]*b[k][j]
            return C

        except: #por si la matriz b es un vector nx1, ya que len(float) es error
            cdef list C = [0]*len(a)
            cdef int i
            for i in range(len(a)):
                cdef int k
                for k in range(len(a[0])):
                    C[i]+=a[i][k]*b[k]
            return C
    cpdef list transpuesta(list m):
        return [list(i) for i in zip(*m)]
    cpdef list menor(list m,int i,int j): #cuando uno tapa la fila i esima con la columna j esima
        return [fila[:j] + fila[j+1:] for fila in (m[:i]+m[i+1:])]
    cpdef float determinante(list m):
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]
        determinant = 0
        cdef int c
        for c in range(len(m)):
            determinant += ((-1)**c)*m[0][c]*determinante(menor(m,0,c))
        return determinant
    cpdef list inversa(list m):
        cdef int determinant = determinante(m)
        if len(m) == 2:
            return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]
        #matriz de cofactores
        cdef list cofactores = []
        cdef int r
        for r in range(len(m)):
            cdef list filacofactores = []
            cdef int c
            for c in range(len(m)):
                cdef list menorr = menor(m,r,c)
                filacofactores.append(((-1)**(r+c)) * determinante(menorr))
            cofactores.append(filacofactores)
        cofactores = transpuesta(cofactores)
        cdef int r
        for r in range(len(cofactores)):
            cdef int c
            for c in range(len(cofactores)):
                cofactores[r][c] = cofactores[r][c]/determinant
        return cofactores
    cpdef list minimos_cuadrados(list t,list y):
        cdef list Xt=[[],[],[]]
        cdef int i
        for i in range(len(t)):
            Xt[0].append(1)
            Xt[1].append(t[i])
            Xt[2].append(t[i]**2)
        X=[list(i) for i in zip(*Xt)]
        cdef list a = mult(mult(inversa(mult(Xt,X)),Xt),y)    
        return a
    cdef list tiempos=list(range(1,len(lc)+1))
    cdef list coeficientes=minimos_cuadrados(tiempos,lc) #y(t)=c+bt+0.5at2 lc son la lista de las posiciones
    cdef float gravedad=coeficientes[-1]*2.0
    cpdef float  SI(float g,float Uespacio,float freq): #Uespacio es el tamano de cada pixel
        return g*(freq**2)*Uespacio #para que de en mm/s^2 si las entradas estan en mm y s
    return SI(gravedad,espacio,frecuencia)
