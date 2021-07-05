import pulp as pl
import random
import copy
import time

#función que leer el archivo y lo guarda en una lista los elementos separador por un caracter de separador
def leer(directorio,separador):
    f=open(directorio,'r')
    C=[]
    for linea in f.readlines():
        linea_v=linea.split(separador)
        if len(linea_v)>1:
            for j in range(len(linea_v)):
                linea_v[j]=float(linea_v[j]) # se usan numero continuos por seguridad 
            C.append(linea_v)
    f.close()
    return C

#recibe un vector de coordenadas y devuelve una matriz de distancias euclidianas
def matriz(coor):
    D=[]
    n=len(coor)
    fila=[]
    for i in range(n):
        for j in range(n):
            fila.append(dist(coor[i],coor[j]))
        D.append(fila)
        fila=[]
    return D

#obtiene la distancia euclidiana entre 2 nodos, a traves de sus coordenadas
def dist(nodo1,nodo2):
    suma1=pow(nodo1[0]-nodo2[0],2)
    suma2=pow(nodo1[1]-nodo2[1],2)
    raiz=round(pow(suma1+suma2,0.5),2)
    return raiz

#Modelo matematico mTSPD
def mTSPD(nombre,tau_T,tau_D,arcos,m,n):# n=nodos #m=camiones
    #arcos=conjunto de rutas posibles de los camiones
    # PARAMETROS
    M=1000000000000
    # CONJUNTOS
    C=[i+1 for i in range(n-1)] #conjunto de clientes
    C0=[i+1 for i in range(n-1)] #conjunto de nodos de salida
    C0.insert(0,'s')
    Cmas=[i+1 for i in range(n-1)] #conjunto de nodos de llegada
    Cmas.append('r')
    CC=[i+1 for i in range(n-1)] #conjunto de todos los nodos
    CC.append('r')
    CC.insert(0,'s')
    F=[(i,j,k) for i in C0 for j in C for k in Cmas] #conjunto de todas las rutas posibles de los drones
    # VARIABLES
    X = pl.LpVariable.dicts("X",arcos,cat='Binary') #variable de desicion de la ruta de los camiones
    Y = pl.LpVariable.dicts("Y",F,cat='Binary') #variable de desicion de la ruta de los drones
    Tl = pl.LpVariable.dicts('Tl',CC,cat='Continuous') #variable de desicion de los tiempos de llegada de los camiones
    Dl = pl.LpVariable.dicts('Dl',CC,cat='Continuous') #variable de desicion de los tiempos de llegada de los drones
    u = pl.LpVariable.dicts('u',C,cat='Continuous') #variable para restringir los subtours
    # FUNCION OBJETIVO
    mdl = pl.LpProblem(nombre, pl.LpMinimize) # modelo de minimizacion
    # 1
    mdl += Tl['r'] #funcion objetivo
    # RESTRICCIONES
    # 2
    for j in C:
        mdl += pl.lpSum(X[(i,j)] for i in C0 if i!=j) + pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if k!=j) for i in C0 if i!=j) == 1
    # 3
    mdl += pl.lpSum(X[('s',j)] for j in Cmas) == m
    # 4
    mdl += pl.lpSum(X[(j,'r')] for j in C0) == m
    # 5
    for i in C:
        for j in C:
            if i!=j:
                mdl += u[i] - u[j] + n*X[(i,j)] + (n-2)*X[(i,j)] <= n-1
    # 6
    for i in C:
        mdl += 1 + (n-2)*X[(i,'r')] + pl.lpSum(X[(j,i)] for j in C if i!=j) <= u[i]
        mdl += u[i] <= n - (n-2)*X[('s',i)] - pl.lpSum(X[(i,j)] for j in C if i!=j)
    # 7
    for j in C:
        mdl += pl.lpSum(X[(i,j)] for i in C0 if i!=j) == pl.lpSum(X[(j,k)] for k in Cmas if k!=j)
    # 8
    for i in C0:
        mdl += pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if i!=k) for j in C if i!=j) <= 1
    # 9
    for k in Cmas:
        mdl += pl.lpSum(pl.lpSum(Y[(i,j,k)] for i in C0 if i!=k) for j in C if j!=k) <= 1
    # 10
    for i in C:
        for j in C:
            for k in Cmas:
                mdl += 2*Y[(i,j,k)] <= pl.lpSum(X[(h,i)] for h in C0 if h!=i) + pl.lpSum(X[(l,k)] for l in C if l!=k)
    # 11
    for j in C:
        for k in Cmas:
            Y[('s',j,k)] <= pl.lpSum(X[(h,k)] for h in C0 if h!=k)
    # 12
    for j in C:
        mdl += pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if k!=j) for i in C0 if i!=j) <= 1 - pl.lpSum(pl.lpSum(Y[(j,a,b)] for b in Cmas if b!=j) for a in C if a!=j)
    # 13
        mdl += pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if k!=j) for i in C0 if i!=j) <= 1 - pl.lpSum(pl.lpSum(Y[(a,b,j)] for a in C0 if a!=j) for b in C if b!=j)
    # 14
    for i in C:
        mdl += Dl[i] >= Tl[i] - M*(1 - pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if i!=k) for j in C if i!=j))
    # 15
        mdl += Dl[i] <= Tl[i] + M*(1 - pl.lpSum(pl.lpSum(Y[(i,j,k)] for k in Cmas if i!=k) for j in C if i!=j))
    # 16
    for k in Cmas:
        mdl += Dl[k] >= Tl[k] - M*(1 - pl.lpSum(pl.lpSum(Y[(i,j,k)] for i in C0 if i!=k) for j in C if j!=k))
    # 17
        mdl += Dl[k] <= Tl[k] + M*(1 - pl.lpSum(pl.lpSum(Y[(i,j,k)] for i in C0 if i!=k) for j in C if j!=k))
    # 18
    for k in Cmas:
        for h in C0:
            mdl += Tl[k] >= Tl[h] + tau_T[(h,k)] - M*(1-X[(h,k)])
    # 19
    for k in Cmas:
        for i in C0:
            for j in C:
                mdl += Dl[k] >= Dl[i] + tau_D[(i,j)] + tau_D[(j,k)] - M*(1-Y[(i,j,k)])
    # 20
    mdl += Dl['s']==0
    # 21
    mdl += Tl['s']==0
    # 22
    for i in Cmas:
        mdl += Dl[i] >= 0
    # 23
        mdl += Tl[i] >= 0
    # 24 
    for i in C:
        for j in C:
            mdl += Y[(i,j,i)] == 0
    # 25
    for i in C:
        mdl += X[(i,i)] == 0    
    # 26
    mdl += X[('s','r')] == 0
    # SOLUCION
    status = mdl.solve() # resolver el problema
    salida=[]
    for v in mdl.variables(): 
        if v.varValue !=0:
            k=b_indices(v) 
            if len(k)>1: # se obtienen los indices de las rutas
                salida.append(k)
    rutas_T,rutas_D=orden(salida,m,n) #se obtienen las rutas de los camiones y drones
    tiempos_T,tiempos_D=tiempo_lleg(m,tau_T,rutas_T,tau_D,rutas_D) #se obtienen los tiempos de llegada de los drones y camiones a los nodos
    T=[]
    T.append(rutas_T)
    T.append(tiempos_T)
    D=[]
    D.append(rutas_D)
    D.append(tiempos_D)
    return T,D,objetivo(tiempos_T,tiempos_D)

#Modelo matematico TSP
def TSP(nombre,tau_T,arcos,n):# n=nodos #m=camiones 
    # PARAMETROS
    M=1000000000000
    # CONJUNTOS
    C=[i+1 for i in range(n-1)] #conjunto de clientes
    C0=[i+1 for i in range(n-1)] #conjunto de nodos de salida
    C0.insert(0,'s')
    Cmas=[i+1 for i in range(n-1)] #conjunto de nodos de llegada
    Cmas.append('r')
    CC=[i+1 for i in range(n-1)] #conjunto de todos los nodos
    CC.append('r')
    CC.insert(0,'s')
    # VARIABLES
    X = pl.LpVariable.dicts("X",arcos,cat='Binary') #variable de desicion de la ruta de los camiones
    Tl = pl.LpVariable.dicts('Tl',CC,cat='Continuous') #variable de desicion de los tiempos de llegada de los camiones
    u = pl.LpVariable.dicts('u',C,cat='Continuous') #variable para evitar los subtours
    # FUNCION OBJETIVO
    mdl = pl.LpProblem(nombre, pl.LpMinimize) # modelo de minimizacion
    # 1
    mdl += Tl['r'] #funcion objetivo
    # RESTRICCIONES
    # en cada nodo llega un solo camion 
    for j in Cmas:
        mdl += pl.lpSum(X[(i,j)] for i in C0 if i!=j) == 1
    # en cada nodo sale un solo camion
    for j in C0:
        mdl += 1 == pl.lpSum(X[(j,k)] for k in Cmas if j!=k)
    # evitar los subtours
    for i in C:
        for j in C:
            if i!=j:
                mdl += u[i] - u[j] + n*X[(i,j)] <= n-1
    # tiempo de llegada de los camiones
    for k in Cmas:
        for h in C0:
            mdl += Tl[k] >= Tl[h] + tau_T[(h,k)] - M*(1-X[(h,k)])
    # el tiempo llegada al inicio es cero
    mdl += Tl['s']==0
    # los tiempos de llegada son positivos
    for i in CC:
        mdl += Tl[i] >= 0
    # no se puede salir y llegar al mismo nodo
    for i in C:
        mdl += X[(i,i)] == 0    
    mdl += X[('s','r')] == 0
    # SOLUCION
    status = mdl.solve() # resolver el problema
    salida=[]
    for v in mdl.variables(): 
        if v.varValue !=0:
            k=b_indices(v)
            if len(k)>1: # se obtienen los indices de las rutas
                salida.append(k)
    ruta_T,D=orden(salida,1,n)
    tiempos_T=tiempo_lleg(1,tau_T,ruta_T)
    rutas_T=[]
    rutas_T.append(ruta_T)
    rutas_T.append(tiempos_T)
    return rutas_T,tiempos_T[0][-1]

#funcion que ordena las rutas de los camiones y drones
#entrada= [1,4],[4,3],[3,5] #salida= [1,4,3,5]
def orden(lista,m,c):
    rutas_T=[]
    rutas_D=[]
    ruta_T=[]
    ruta_D=[]
    n=len(lista)
    indices=[]
    clientes=[i+1 for i in range(c)]
    clientes.insert(0,'s')
    clientes.append('r')
    for i in range(n):
        if len(lista[i])==2:
            rutas_T.append(lista[i])
        else:
            rutas_D.append(lista[i])
    indices.append(buscar(rutas_T,'s',m))
    for i in clientes[1:c+1]:
        indices.append(buscar(rutas_T,i,2))
    indices.append(buscar(rutas_T,'r',m))
    vis=[]
    for i in range(m):
        vis.append('s')
        sig=rutas_T[i][1]
        vis.append(sig)
        while sig!='r':
            sig=rutas_T[indices[sig][0]][1]
            vis.append(sig)
        ruta_T.append(vis)
        vis=[]
    return ruta_T,rutas_D

#funcion que busca el elemento k en la lista rutas y devuelve lim indices(1ra dimension) donde esta
def buscar(rutas,k,lim):
    indices=[]
    for i in range(len(rutas)):
        if k==rutas[i][0]:
            indices.append(i)
        if len(indices)==lim:
            break
    return indices

#funcion que busca el elemento k en la lista rutas y duevuelve el ultimo indice, si det=True devuelve el otro indice (2da dimension)
def buscar2(rutas,k,det=True):
    h=None
    j=None
    for i in range(len(rutas)):
        if k in rutas[i]:
            if det==True:
                try:
                    j=rutas[i].index(k)
                except:
                    j=None
            h=i
            break
    if det==True:
        t=h,j
    else:
        t=h
    return t

#Funcion que busca el elemento k en la lista rutas y devuelve el ultimo n par de indices encontrados
def buscar3(rutas,k,n):
    h=None
    j=None
    det=True
    cuenta=0
    for i in range(len(rutas)):
        if k in rutas[i]:
            try:
                j=rutas[i].index(k)
            except:
                j=None
            h=i
            cuenta+=1
            if cuenta==n:
                break
    t=h,j
    return t

#funcion que separa una lista en 2, eliminando n elementos aleatorios de X transfiriendolos a lista
def elim(X,n):
    lista=[]
    for i in range(n):
        indice=random.randint(0,len(X)-1)
        lista.append(X.pop(indice))
    return X,lista

#funcion que crea un cromosoma para el GA con n=nodos y m=camiones
def cromosoma(n,m):
    nodos=[i+1 for i in range(n-1)]
    rutas=[]
    tam=len(nodos)
    i=0
    while i<m-1 and tam!=1:
        l_max=tam-(m-1-i)
        l_ruta=random.randint(1,l_max)
        nodos,ruta_i=elim(nodos,l_ruta)
        ruta_i.insert(0,'s')
        ruta_i.append('r')
        rutas.append(ruta_i)
        tam=len(nodos)
        i+=1
    ruta_i=nodos
    ruta_i.insert(0,'s')
    ruta_i.append('r')
    rutas.append(ruta_i)
    return rutas

#obtiene el valor aptitud de un cromosoma
def aptitud(tau,ruta,m):
    t_max=0
    for i in range(m):
        tiempo=0
        for j in range(len(ruta[i])-1):
            tiempo=tiempo+tau[ruta[i][j],ruta[i][j+1]]
        if tiempo>t_max:
            t_max=tiempo
    t_max=round(t_max,2)
    return t_max

#obtiene los indices de la variable v
def b_indices(v):
    name=v.name
    name=name.strip('X_(')
    name=name.strip('Y_(')
    name=name.strip(')')
    k=name.split(',_')
    for i in range(len(k)):
        try:
            k[i]=int(k[i])
        except ValueError:
            k[i]=str(k[i])
            k[i]=k[i].strip("'")
    return k

#busca el minimo de lista, sin consideran los indices de lista_indices
def min_sin(lista,lista_indices,M):
    j=None
    mini=M
    for i in range(len(lista)):
        if i not in lista_indices and lista[i]<mini:
            j=i
            mini=lista[i]
    return j

#busca los n elementos minimos de lista
def mintop(lista,n,M):
    min=[]
    for i in range(n):
        indice=min_sin(lista,min,M)
        min.append(indice)
    return min

#conjunto de funciones que cruzan 2 cromosomas y que cumplan con las restricciones del mTSP
def crucein(lista1,lista2):
    aux1=[]
    aux2=[]
    n=len(lista1)
    ref=[i+1 for i in range(n)]
    for i in range(0,n,2):
        if lista1[i] not in aux1:
            aux1.append(lista1[i])
        if i!=n-1 and lista2[i+1] not in aux1:
            aux1.append(lista2[i+1])
    for i in range(0,n,2):
        if lista2[i] not in aux2:
            aux2.append(lista2[i])
        if i!=n-1 and lista1[i+1] not in aux2:
            aux2.append(lista1[i+1])
    for i in range(n):
        if ref[i] not in aux1:
            aux1.append(ref[i])
        if ref[i] not in aux2:
            aux2.append(ref[i])
    return aux1,aux2

def cruce(papa,mama,m):
    aux1=[]
    aux2=[]
    nn1=[0]
    nn2=[0]
    for i in range(m):
        n1=len(papa[i])
        nn1.append(n1-2+nn1[i])
        n2=len(mama[i])
        nn2.append(n2-2+nn2[i])
        for j in range(1,n1-1):
            aux1.append(papa[i][j])
        for j in range(1,n2-1):
            aux2.append(mama[i][j])
    aux1,aux2=crucein(aux1,aux2)
    aux1=separar(aux1,nn1)
    aux2=separar(aux2,nn2)
    return aux1,aux2

def separar(lista,separador):
    salida=[]
    for i in range(len(separador)-1):
        aux=lista[separador[i]:separador[i+1]]
        aux.insert(0,'s')
        aux.append('r')
        salida.append(aux)
    return salida

#conjunto de funciones que mutan a los cromosomas
def mutacion1(c,m):
    cromosoma=copy.deepcopy(c)
    for i in range(m):
        n=len(cromosoma[i])
        if n-2>1:
            indice1=random.randint(1,n-2)
            indice2=random.randint(1,n-2)
            while indice2==indice1:
                indice2=random.randint(1,n-2)
            aux=cromosoma[i][indice1]
            cromosoma[i][indice1]=cromosoma[i][indice2]
            cromosoma[i][indice2]=aux
    return cromosoma

def mutacion2(c,m):
    cromosoma=copy.deepcopy(c)
    if m>1:
        r1=random.randint(0,m-1)
        r2=random.randint(0,m-1)
        while r2==r1:
            r2=random.randint(0,m-1)
        n1=len(cromosoma[r1])
        n2=len(cromosoma[r2])
        indice1=random.randint(1,n1-2)
        indice2=random.randint(1,n2-2)
        aux=cromosoma[r1][indice1]
        cromosoma[r1][indice1]=cromosoma[r2][indice2]
        cromosoma[r2][indice2]=aux
    return cromosoma

#heuristica genetica que resuelve el problema como un mTSP
def GA(tau_T,m,n,NumGen,PopSize,pc,pm,ps):
    M=1000000000
    indices=[i for i in range(PopSize)]
    # crear poblacion inicial
    Pob=[]
    sol=[]
    bs_Z=M
    for i in range(PopSize):
        cro=cromosoma(n,m)
        Z=aptitud(tau_T,cro,m)
        Pob.append(cro)
        sol.append(Z)
        #encontrar la mejor solucion dentro de la poblacion inicial
        if Z<bs_Z:
            bs_Z=Z
            bs=copy.deepcopy(cro)
    for i in range(NumGen):
        #generacion de nuevos individuos
        for j in range(PopSize-len(Pob)):
            cro=cromosoma(n,m)
            Z=aptitud(tau_T,cro,m)
            Pob.append(cro)
            sol.append(Z)
        #Cruzamiento
        n_sel=int(pc*PopSize)
        i_sel=random.sample(indices,n_sel)
        for i in range(n_sel-1):
            [hijo1,hijo2]=cruce(Pob[i_sel[i]],Pob[i_sel[i+1]],m)
            Pob.append(hijo1)
            sol.append(aptitud(tau_T,hijo1,m))
            Pob.append(hijo2)
            sol.append(aptitud(tau_T,hijo2,m))
        # Mutaciones
        n_mut=int(pm*PopSize)
        i_mut=random.sample(indices,n_mut)
        for i in range(n_mut):
            mut=mutacion1(Pob[i],m)
            mut2=mutacion2(Pob[i],m)
            Pob.append(mut)
            Pob.append(mut2)
            sol.append(aptitud(tau_T,mut,m))
            sol.append(aptitud(tau_T,mut2,m))
        #buscar al mejor
        for i in range(len(Pob)):
            if sol[i]<bs_Z:
                bs_Z=sol[i]
                bs=copy.deepcopy(Pob[i])
        #mantener al ps% mejor de la poblacion para la siguiente generacion
        indices_pob=mintop(sol,int(ps*PopSize),M)
        Pob=[Pob[i] for i in indices_pob]
        sol=[sol[i] for i in indices_pob]
    return bs

#funcion que obtiene la matriz de tiempos para solamente los camiones
def obtener_tau(dist,vel_T):
    vel_D=vel_T*1.5
    tau_T=[]
    arcos=[]
    tau_D=[]
    for i in range(len(dist)):
        for j in range(len(dist[0])):
            if i==0 and j==0:
                tau_T.append(dist[i][j]/vel_T)
                tau_D.append(dist[i][j]/vel_T)
                arcos.append(('s','r'))
            elif i==0:
                tau_T.append(dist[i][j]/vel_T)
                tau_D.append(dist[i][j]/vel_T)
                arcos.append(('s',j))
            elif j==0:
                tau_T.append(dist[i][j]/vel_T)
                tau_D.append(dist[i][j]/vel_T)
                arcos.append((i,'r'))
            else:
                tau_T.append(dist[i][j]/vel_T)
                tau_D.append(dist[i][j]/vel_T)
                arcos.append((i,j))
    tau_T=dict(zip(arcos,tau_T))
    tau_D=dict(zip(arcos,tau_D))
    return tau_T,tau_D,arcos

# Funcion que obtiene la matriz de tiempos para los camiones y los drones
def obtener_tau2(dist,vel_T):
    tau_T=[]
    arcos=[]
    for i in range(len(dist)):
        for j in range(len(dist[0])):
            if i==0 and j==0:
                tau_T.append(dist[i][j]/vel_T)
                arcos.append(('s','r'))
            elif i==0:
                tau_T.append(dist[i][j]/vel_T)
                arcos.append(('s',j))
            elif j==0:
                tau_T.append(dist[i][j]/vel_T)
                arcos.append((i,'r'))
            else:
                tau_T.append(dist[i][j]/vel_T)
                arcos.append((i,j))
    tau_T=dict(zip(arcos,tau_T))
    return tau_T,arcos

#funcion que obtiene los tiempos de llegada de los camiones y drones considerando que los tiempos de llegada de los camiones y drones son los mismos para el mismo nodo
def tiempo_lleg(m,tau_T,ruta_T,tau_D=None,ruta_D=None):
    tiempos=[]
    tiempos_T=[]
    tiempos_D=[]
    for i in range(m):
        tiempo=[0]
        for j in range(1,len(ruta_T[i])):
            tiempo.append(round(tau_T[ruta_T[i][j-1],ruta_T[i][j]]+tiempo[j-1],2))
        tiempos_T.append(tiempo)
    if ruta_D!=None:
        for i in range(len(ruta_D)):
            h,k=buscar2(ruta_T,ruta_D[i][0],True)
            tiempo=[tiempos_T[h][k]]
            for j in range(1,len(ruta_D[i])):
                tiempo.append(round(tau_D[ruta_D[i][j-1],ruta_D[i][j]]+tiempo[j-1],2))
            tiempos_D.append(tiempo)
        for i in range(len(ruta_D)):
            if ruta_D[i][-1]!='r':
                h,k=buscar2(ruta_T,ruta_D[i][-1],True)
                if tiempos_D[i][-1]>tiempos_T[h][k]:
                    tiempos_T[h][k]=tiempos_D[i][-1]
                    tiempos_T[h]=recal_tiempos(ruta_T[h],tiempos_T[h],k,tau_T)
                else:
                    tiempos_D[i][-1]=tiempos_T[h][k]
        tiempos.append(tiempos_T)
        tiempos.append(tiempos_D)
    else:
        tiempos=tiempos_T
    return tiempos

#funcion que elimina un nodo de una ruta de camiones y/o drones
def node_removal(T_R,D_P,j,tau_T,tau_D,camiones,Lch_N,La_N):
    h1,k1=buscar2(T_R[0],j,True)
    #camiones
    if h1!=None and k1!=None:
        camion=True
        T_R[0][h1].pop(k1)
        T_R[1]=tiempo_lleg(camiones,tau_T,T_R[0])
        t1=h1,k1
        if len(D_P)>0:
            h2,k2=buscar2(D_P[0],j,True)
            if h2!=None and k2!=None and k2!=1:
                i=2
                while h2!=None and k2!=None and k2!=1:
                    if k2==0:
                        if T_R[0][h1][k1-1] not in Lch_N:
                            aux=T_R[0][h1][k1-1]
                        else:
                            aux='s'
                        D_P[0][h2][k2]=copy.deepcopy(aux)
                        Lch_N.append(aux)
                    elif k2==2:
                        if T_R[0][h1][k1] not in La_N:
                            aux=T_R[0][h1][k1]
                        else:
                            aux='r'
                        D_P[0][h2][k2]=copy.deepcopy(aux)
                        La_N.append(aux)
                    t2=h2,k2
                    t=[t1,t2]
                    h2,k2=buscar3(D_P[0],j,i)
                    i+=1
                T_R[1],D_P[1]=tiempo_lleg(camiones,tau_T,T_R[0],tau_D,D_P[0])
            else:
                t=[t1,(None,None)]
        else:
            t=[t1,(None,None)]
            D_P=[[],[]]
    return T_R,D_P,t

#funcion que calcula el numero de drones y los nodos donde salen y aterrizan
def cal_drones(D_P):
    if len(D_P)!=0:
        NCS=len(D_P[0])
    else:
        NCS=0
    Lch_N=[]
    La_N=[]
    for i in range(NCS):
        if D_P[0][i][0]!='s' and D_P[0][i][0] not in Lch_N:
            Lch_N.append(D_P[0][i][0])
        if D_P[0][i][2]!='r' and D_P[0][i][2] not in La_N:
            La_N.append(D_P[0][i][2])
    return NCS,Lch_N,La_N

#funcion que permite coincidir los tiempos de llegada de los drones y camiones al mismo nodo
def recal_tiempos(ruta,tiempos,t1_1,tau_T):
    for i in range(t1_1+1,len(ruta)):
        tiempos[i]=round(tau_T[ruta[i-1],ruta[i]]+tiempos[i-1],2)
    return tiempos

#funcion que agrega un nodo como dron a la misma ruta (DronAddedSameRoute)
def DA_SR(T_R_,D_P_,j,tau_T,tau_D,camiones,t1,t2,Lch_N,La_N):
    T_R=copy.deepcopy(T_R_)
    D_P=copy.deepcopy(D_P_)
    nodo_ant=T_R[0][t1[0]][t1[1]-1]
    i=2
    while nodo_ant in Lch_N:
        nodo_ant=T_R[0][t1[0]][t1[1]-i]
        i+=1
    nodo_sig=T_R[0][t1[0]][t1[1]]
    i=1
    while nodo_sig in La_N:
        nodo_sig=T_R[0][t1[0]][t1[1]+i]
        i+=1
    ruta_dron=[nodo_ant,j,nodo_sig]
    D_P[0].append(ruta_dron)
    tpo=tiempo_lleg(camiones,tau_T,T_R[0],tau_D,D_P[0]) # tpo tiene los tiempos de los camiones[0] y drones[1]
    T_R[1]=tpo[0]
    D_P[1]=tpo[1]
    return T_R,D_P

#funcion que encuentra el nodo mas cercano a j y que no sea de la ruta t1
def mas_cercano(T_R,tau_D,j,t1,camiones):
    min=10000000000000
    for i in range(camiones):
        if i!=t1[0]:
            for k in range(len(T_R[0][i])-1):
                if tau_D[j,T_R[0][i][k+1]]<min:
                    ind_min=i,k+1
                    min=tau_D[j,T_R[0][i][k+1]]
    return ind_min

#funcion que encuentra un numero aleatorio que no sea sin, y otro indice(2da dimension) del vector[i]
def aleatorio_sin(inicio,fin,sin,T_R_0):
    i=random.randint(inicio,fin-1)
    while i==sin:
        i=random.randint(inicio,fin-1)
    j=random.randint(1,len(T_R_0[i])-1)
    return i,j

#funcion que calcula el valor objetivo de todas las rutas
def objetivo(T_T,T_D):
    lista=[]
    for i in range(len(T_T)):
        lista.append(T_T[i][-1])
    for i in range(len(T_D)):
        lista.append(T_D[i][-1])
    return max(lista)

#funcion que agrega un nodo como dron a una ruta distinta (DronAddedDiffRoute)
def DA_DR(T_R_,D_P_,j,tau_T,tau_D,camiones,t1,t2,Lch_N,La_N):
    T_R=copy.deepcopy(T_R_)
    D_P=copy.deepcopy(D_P_)
    w,q=mas_cercano(T_R,tau_D,j,t1,camiones)
    nodo_ant=T_R[0][t1[0]][t1[1]-1]
    i=2
    while nodo_ant in Lch_N:
        nodo_ant=T_R[0][t1[0]][t1[1]-i]
        i+=1
    nodo_sig=T_R[0][w][q]
    i=1
    while nodo_sig in La_N:
        nodo_sig=T_R[0][w][q+i]
        i+=1
    ruta_dron=[nodo_ant,j,nodo_sig]
    D_P[0].append(ruta_dron)
    tpo=tiempo_lleg(camiones,tau_T,T_R[0],tau_D,D_P[0]) # tpo tiene los tiempos de los camiones[0] y drones[1]
    T_R[1]=tpo[0]
    D_P[1]=tpo[1]
    return T_R,D_P

#funcion que agrega un nodo a una ruta distinta (TruckInsertion)
def TI(T_R_,D_P_,j,tau_T,tau_D,camiones,t1,t2):
    T_R=copy.deepcopy(T_R_)
    D_P=copy.deepcopy(D_P_)
    w,q=mas_cercano(T_R,tau_D,j,t1,camiones)
    T_R[0][w].insert(q,j)
    tpo=tiempo_lleg(camiones,tau_T,T_R[0],tau_D,D_P[0]) # tpo tiene los tiempos de los camiones[0] y drones[1]
    T_R[1]=tpo[0]
    D_P[1]=tpo[1]
    return T_R,D_P

#heuristica mixta que resuelve el problema como mTSPD
def ADI(tau_T,tau_D,camiones,clientes):
    M=10000000000000
    #1ra Fase
    tours=GA(tau_T,camiones,clientes,50,50,0.5,0.5,0.5)
    t_lleg=tiempo_lleg(camiones,tau_T,tours)
    #inicializacion de parametros
    max_NCS=0
    for i in range(camiones):
        max_NCS=max_NCS+int((len(tours[i])-2)/2)
    Lch_N=[]
    La_N=[]
    T_R=[tours,t_lleg]
    ahorro=1
    cand_list=[i+1 for i in range(clientes-1)]
    D_P=[]
    delta_T_S=0
    T_D=objetivo(T_R[1],D_P)
    T_R_best=T_R
    D_P_best=D_P
    Z_best=objetivo(T_R[1],D_P)
    while ahorro >0:
        T_R_dic={}
        D_P_dic={}
        deltas={}
        T_R=copy.deepcopy(T_R_best)
        D_P=copy.deepcopy(D_P_best)
        for j in cand_list:
            h=buscar2(T_R[0],j,False)
            if h != None:
                select_tour=T_R[0][h]
            else:
                h=buscar2(D_P[0],j,False)
                select_tour=D_P[0][h]
            if len(select_tour[1:len(select_tour)-1]) < 2:
                delta_T_S=-M
            else:
                NCS,Lch_N,La_N=cal_drones(D_P)
                T_R,D_P,[t1,t2]=node_removal(T_R,D_P,j,tau_T,tau_D,camiones,Lch_N,La_N)
                NCS,Lch_N,La_N=cal_drones(D_P)
                T=[]
                if NCS < max_NCS:
                    T_R1,D_P1=DA_SR(T_R,D_P,j,tau_T,tau_D,camiones,t1,t2,Lch_N,La_N)
                    tm1=objetivo(T_R1[1],D_P1[1])
                    T.append(tm1)
                    T_R2,D_P2=DA_DR(T_R,D_P,j,tau_T,tau_D,camiones,t1,t2,Lch_N,La_N)
                    tm2=objetivo(T_R2[1],D_P2[1])
                    T.append(tm2)
                else:
                    tm12=M
                    T.append(tm12)
                    T.append(tm12)
                T_R3,D_P3=TI(T_R,D_P,j,tau_T,tau_D,camiones,t1,t2)
                tm3=objetivo(T_R3[1],D_P3[1])
                T.append(tm3)
                ###
                ind_T=T.index(min(T))
                T_I=T[ind_T]
                delta_T_S=T_D-T_I
                if ind_T==0:
                    T_R=T_R1
                    D_P=D_P1
                elif ind_T==1:
                    T_R=T_R2
                    D_P=D_P2
                elif ind_T==2:
                    T_R=T_R3
                    D_P=D_P3
            T_R_dic[j]=copy.deepcopy(T_R)
            D_P_dic[j]=copy.deepcopy(D_P)
            deltas[j]=copy.deepcopy(delta_T_S)
        lista_objetivo=[deltas[i] for i in cand_list]
        ahorro=max(lista_objetivo)
        indice=lista_objetivo.index(ahorro)
        Z_ahorro=objetivo(T_R_dic[cand_list[indice]][1],D_P_dic[cand_list[indice]][1])
        if ahorro>0:
            T_R_best=T_R_dic[cand_list[indice]]
            D_P_best=D_P_dic[cand_list[indice]]
            if len(D_P_best)>0:
                Z_best=objetivo(T_R_best[1],D_P_best[1])
                T_D=objetivo(T_R[1],D_P[1])
            else:
                Z_best=objetivo(T_R_best[1],[])
                T_D=objetivo(T_R[1],[])
        elim=cand_list.pop(indice)
        # print('---------------------------------------------------------')
        # print(cand_list)
        # print(elim)
        # print(T_R_best)
        # print(D_P_best)
        # print(Z_best)
    return T_R_best,D_P_best,Z_best

#DIRECCION ABSOLUTA DE LA CARPETA "instancias"
directorio="C:\\Users\\usuario\\Documents\\UNIVERSIDAD\\SEXTO\\TALLER DE OPTIMIZACION\\TRABAJO\\instancias\\"

#mTSPD vs TSP
r1=open(directorio+'ciudades\\resultados_ciudades.txt','w')
r1_=open(directorio+'ciudades\\resultados_ciudades_rutas.txt','w')
r1.write('Problema,Camiones,Objetivo mTSPD,Tiempo de Ejecución mTSPD,Objetivo TSP,Tiempo de Ejecución TSP\n')
for i in ['calama','concepcion','santiago']:
    dist_T=leer(directorio+'ciudades\\'+i+'_camion.txt',',')
    dist_D=leer(directorio+'ciudades\\'+i+'_dron.txt',',')
    clientes=len(dist_T)
    vel_T=10 #10 metros por segundo
    tau_T,arcos=obtener_tau2(dist_T,vel_T)
    tau_D,arcos=obtener_tau2(dist_D,vel_T*1.5)
    inicio2=time.time()
    rutas_T_TSP,obj_TSP=TSP(i,tau_T,arcos,clientes)
    fin2=time.time()
    t_eje_TSP=round(fin2-inicio2,4)
    r1_.write('-----------------------------------------------------------------------\n')
    r1_.write('Problema'+i+':\n')
    r1_.write('TSP= '+str(rutas_T_TSP[0])+'\nTpo= '+str(rutas_T_TSP[1])+'\n')
    r1_.write('mTSPD:')
    for camiones in range(2,4):
        inicio2=time.time()
        rutas_T_mTSPD,rutas_D_mTSPD,obj_mTSPD=mTSPD(i,tau_T,tau_D,arcos,camiones,clientes)
        fin2=time.time()
        t_eje_mTSPD=round(fin2-inicio2,4)
        r1_.write('Camiones:'+str(camiones)+'\n')
        r1_.write('Camiones:'+str(rutas_T_mTSPD[0])+'\nTpo= '+str(rutas_T_mTSPD[1])+'\n')
        r1_.write('Drones:'+str(rutas_D_mTSPD[0])+'\nTpo= '+str(rutas_D_mTSPD[1])+'\n')
        r1_.write('-----------------------------------------------------------------------\n')
        r1.write(str(i)+','+str(camiones)+','+str(obj_mTSPD)+','+str(t_eje_mTSPD)+','+str(obj_TSP)+','+str(t_eje_TSP)+'\n')
r1.close()
r1_.close()
#HEURISTICA - TSPLIB
f=open(directorio+'tsplib\\nombres_tsp_lib.txt','r')
r2=open(directorio+'tsplib\\resultados_tsp_lib.txt','w')
r2.write('Problema,Camiones,Objetivo Promedio,Mejor Objetivo,Tiempo de Ejecución Promedio\n')
for i in f.readlines():
    nombre=i.strip('\n')
    coor=leer(directorio+'tsplib\\'+nombre,',')
    ma_dist=matriz(coor)
    clientes=len(ma_dist)
    vel_T=10 # 10 metros por segundo
    tau_T,tau_D,arcos=obtener_tau(ma_dist,vel_T)
    for camiones in range(2,6):
        rutas_camiones=[]
        rutas_drones=[]
        tiempo_final=[]
        tiempo_ejecucion=[]
        for h in range(5):
            inicio=time.time()
            rc,rd,tf=ADI(tau_T,tau_D,camiones,clientes)
            fin=time.time()
            te=fin-inicio
            rutas_camiones.append(rc)
            rutas_drones.append(rd)
            tiempo_final.append(tf)
            tiempo_ejecucion.append(te)
        #mejor solucion
        mejor_tiempo_final=min(tiempo_final)
        ind_mejor=tiempo_final.index(mejor_tiempo_final)
        mejor_tiempo_ejecucion=tiempo_ejecucion[ind_mejor]
        mejor_rutas_camiones=rutas_camiones[ind_mejor]
        mejor_rutas_drones=rutas_drones[ind_mejor]
        #promedio
        prom_tiempo_final=round(sum(tiempo_final)/len(tiempo_final),4)
        prom_tiempo_ejecucion=round(sum(tiempo_ejecucion)/len(tiempo_ejecucion),4)
        #escritura
        r2.write(nombre+','+str(camiones)+','+str(prom_tiempo_final)+','+str(mejor_tiempo_final)+','+str(prom_tiempo_ejecucion)+'\n')
f.close()
r2.close()
