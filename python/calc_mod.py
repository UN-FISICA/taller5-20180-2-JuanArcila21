from numpy import where
import scipy.misc as sp
import scipy.ndimage as nd
from scipy import ndimage
import numpy as np
from numpy.linalg import inv
def calc(img,frecuencia,espacio):
    h,c,w=img.shape #h para filas c para columnas w para numero de matrices 1 con cada color
    ima=np.zeros((h,c)) #matriz con 1 solo color del mismo tamaño que img
    for i in range(h):
        for j in range(c):
            promedio=(int(img[i][j][0])+int(img[i][j][1])+int(img[i][j][2]))/3
            ima[i][j]=int(promedio) #valor de gris igual al promedio de RGB
    sp.imsave("Bolas1BN.png",ima) #guarda la nueva imagen en escala de grises
    imgblancoynegro=sp.imread("Bolas1BN.png") #abre la imagen en escala de grises 
    if imgblancoynegro[0][0]<128:
        imb=where(imgblancoynegro<80,0,255) #binariza la imagen
    else:
        imb=where(imgblancoynegro>200,0,255) 
    imb=nd.median_filter(imb,(2,2)) #la filtra para quitar zonas claras pequenas que no son bolas
    imc,n=nd.label(imb) #el elemento 0 es la matriz, el 1 es el numero de objetos contados
    cm=ndimage.measurements.center_of_mass(imb,imc,list(range(1,n+1))) #halla los centros de masa
    lc=[] #lista coordenadas
    for i in range(len(cm)):
        lc.append(cm[i][0])
    lc.remove(min(lc))
    def minimos_cuadrados(t,y):
        f=[]
        f.append(lambda x:np.ones_like(x))
        f.append(lambda x:x)
        f.append(lambda x:x**2)
        Xt=[]
        
        for fu in f:
            Xt.append(fu(t))
        Xt= np.array(Xt)
        X=Xt.transpose()
        
        a = np.dot(np.dot(inv(np.dot(Xt,X)),Xt),y)
        return a
    tiempos=np.linspace(1,n-1,n-1) #n-1 datos para n bolas porque la primera fue quitada
    posiciones=np.array(lc)
    coeficientes=minimos_cuadrados(tiempos,posiciones) #y(t)=c+bt+0.5at2
    gravedad=coeficientes[-1]*2.0
    #print (gravedad) #en pixeles/unidades de tiempo al cuadrado
    def SI(g,Uespacio,freq): #Uespacio es el tamano de cada pixel
        return g*(freq**2)*Uespacio #para que de en mm/s^2 si las entradas estan en mm y s
    return SI(gravedad,espacio,frecuencia)
