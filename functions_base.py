# email for communication: prioyutipradhan@gmail.com
import networkx as nx                                              
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
from numpy import linalg as LA
import random

def print_evec_entries(A):
  e_val, e_vec = LA.eigh(A)
  d = np.shape(e_vec)   
  n = d[0]  
  Ev = abs(e_vec[:,n-1])
  for i in Ev:
    print '%0.4f' %i 
  del e_val
  del e_vec
  del Ev

def store_undirected_network(A, N, path):
  f_out = open(path,'wb')
  for i in range(N-1):
    for j in range(i+1,N,1):
      if A[i,j] == 1:
        #print i+1,j+1
        f_out.write(str(i+1)+' '+str(j+1)+'\n')
  f_out.close()

def read_network(filename):

   fh = open(filename,"rb")
    
   G = nx.Graph()
   for line in fh:
     x,y = line.split()
     x1 = int(x)
     y1 = int(y) 
     G.add_edge(x1, y1)
   
   #N = len(nx.nodes(G))

   return G

def print_network(A,N):
  for i in range(N):
    for j in range(N):
      print("%d" %A[i,j]),
    print
  print      

def Edge_Rewiring(a, N, Edges):
  T = ((N*(N-1))/2) - Edges
  #print N, Edges, T
  k = 0
  l = 0
  e1 = np.zeros((Edges,2),dtype=int)
  e2 = np.zeros((T,2),dtype=int)
  for i in range(N-1): 
    for j in range(i+1,N,1):
     # print i,j
     if a[i,j] == 1:
        e1[k,0] = i
        e1[k,1] = j
        k = k + 1
     else:
        e2[l,0] = i
        e2[l,1] = j
        l = l + 1
  
  a_New = np.zeros((N,N),dtype=int)
  #print(a)
  while True:
   a_New = np.copy(a)
   #print(a_New)
   r1 = random.randint(0,k-1)
   a_New[e1[r1,0],e1[r1,1]] = 0  
   a_New[e1[r1,1],e1[r1,0]] = 0
    
   r0 = random.randint(0,l-1) 
   a_New[e2[r0,0],e2[r0,1]] = 1 
   a_New[e2[r0,1],e2[r0,0]] = 1
   
   #-------checking for connectedness----------
   G = nx.from_numpy_matrix(a_New)
   #print(nx.is_connected(G))
   #a_New 
   if nx.is_connected(G):  
     del e1 
     del e2
     return a_New 

def main_rewiring(G,Iter):
  Edges = nx.number_of_edges(G)
  print 'Number of edges: ', Edges
  N = len(G.nodes())
  print N
  A = nx.to_numpy_matrix(G, dtype=int)
  print A.shape
  e_val1, e_vec1 = LA.eigh(A)
  print e_vec1[:,N-1].shape
  ev1 = e_vec1[:,N-1]
  ev2 = e_vec1[:,N-2]
  fh = open('dot_product.txt','wb')
  #c1 = sum(np.power(Ev1,4))
  d1 = 0
  d2 = 0
  flag = 0
  fh.write(str(d1)+' '+str(d2)+' '+str(e_val1[N-1])+' '+str(e_val1[N-2])+' '+str(e_val1[N-3])+'\n')
  del e_vec1
  del e_val1
  #del ev1
  #del ev2  

  for iter in range(2,Iter+1,1): 
     flag = 0
     A_New = np.zeros((N,N),dtype = int)
     A_New = Edge_Rewiring(A, N, Edges)
     print A_New.shape
     e_val, e_vec = LA.eigh(A_New)
     #print e_vec[:,N-1].shape
     Ev1 = e_vec[:,N-1]
     Ev2 = e_vec[:,N-2] 
     E1 = Ev1.reshape((N, 1))
     E2 = Ev2.reshape((N, 1))
     #c2 = sum(np.power(Ev,4))
     print E2.shape
     print E1.shape
     d1 = np.dot(ev1,E1)
     d2 = np.dot(ev2,E1)
     
     ev1 = Ev1
     ev2 = Ev2   
      
     '''if c2 > c1:
       A = np.copy(A_New)
       c1 = c2
       flag = 1'''
     #print(c1) 
     fh.write(str(d1)+' '+str(d2)+' '+str(e_val[N-1])+' '+str(e_val[N-2])+' '+str(e_val[N-3])+'\n')      
     del A_New
     del e_vec
     del e_val
     del Ev1
     del Ev2 
  fh.close()
  #G = nx.from_numpy_matrix(A)
  #f_out = open(file_out,'wb')
  #nx.write_edgelist(G, f_out, data=False)
  del A
  #G.clear()
  print 'Completes...'

def get_ipr(A):
  e_val, e_vec = LA.eigh(A)
  d = np.shape(e_vec)   
  n = d[0]   
  lam1 = e_val[n-1]  
  lam2 = e_val[n-2]   
  lam3 = e_val[n-3] 
  Ev1 = e_vec[:,n-1]
  Ev2 = e_vec[:,n-2]
  Ev3 = e_vec[:,0]
          
  IPR1 = 0.0;
  IPR2 = 0.0;
  IPR3 = 0.0;
  for i in range(n):
    IPR1 = IPR1 + pow(Ev1[i],4);
    IPR2 = IPR2 + pow(Ev2[i],4); 
    IPR3 = IPR3 + pow(Ev3[i],4); 
  #print 'IPR of x1, and x2: ',IPR1[0,0], IPR2[0,0],IPR3[0,0]
  del e_val
  del e_vec
  del Ev1
  del Ev2
  return IPR1,IPR2,IPR3,lam1,lam2,lam3
  
def eigval_bound(A):
  deg = sum(A)
  k_max = max(deg)
  k_mean = np.mean(deg);
  k_var = np.var(deg);
  E_lambda = k_max/np.sqrt(k_max-(k_var/k_mean)+1);
  print 'Estimated lambda1: ',E_lambda  
  
def create_connected_ER_random_network(s,p):
   Flag = False
   while not Flag:
     G = nx.fast_gnp_random_graph(s, p)
     Flag = nx.is_connected(G)
     print Flag
   return G

def combining_components(G1,G2,s1,s2,N):  
   M = np.zeros((N,N), dtype=int)
   M1 = np.zeros((s1,s1), dtype=int)
   M2 = np.zeros((s2,s2), dtype=int)

   Edge_list1 = nx.edges(G1)
   for edge in Edge_list1:
     x = edge[0]
     y = edge[1]
     M1[x][y] = 1
     M1[y][x] = 1

   Edge_list2 = nx.edges(G2)
   for edge in Edge_list2:
     x = edge[0]
     y = edge[1]
     M2[x][y] = 1
     M2[y][x] = 1

   print 'Number of edges in random graph: ', nx.number_of_edges(G1)
   ipr1,ipr2,ipr3,lam1,lam2,lam3 = get_ipr(M1)                
   print 's1: %d ipr1: %0.4f ipr2: %0.4f ipr3: %0.4f lam1: %0.4f lam2: %0.4f lam3: %0.4f\n'% (s1,ipr1,ipr2,ipr3,lam1,lam2,lam3)                                                                                                        

   print 'Number of edges in random graph: ', nx.number_of_edges(G2)
   ipr1,ipr2,ipr3,lam1,lam2,lam3 = get_ipr(M2)                
   print 's2: %d ipr1: %0.4f ipr2: %0.4f ipr3: %0.4f lam1: %0.4f lam2: %0.4f lam3: %0.4f\n'% (s2,ipr1,ipr2,ipr3,lam1,lam2,lam3)                                                                                                        

 
   for i in range(s1):
     for j in range(s1):
        M[i][j] = M1[i][j]

   M[s1-1][s1] = 1
   M[s1][s1-1] = 1 
   M[s1][s1+1] = 1
   M[s1+1][s1] = 1

   i = s1 + 1
   j = s1 + 1
   m = 0
   n = 0

   while i<N and m<s2:
      j = s1+1
      n = 0
      while j<N and n<s2:
        M[i][j] = M2[m][n]
        #print i,j,m,n
        j = j + 1
        n = n + 1
      i = i + 1
      m = m + 1

   G = nx.from_numpy_matrix(M)
   print 'Number of edges in random-random graph combined: ', nx.number_of_edges(G)
   ipr1,ipr2,ipr3,lam1,lam2,lam3 = get_ipr(M)                
   print 'N: %d\n ipr1: %0.4f\n ipr2: %0.4f\n ipr3: %0.4f\n lam1: %0.4f\n lam2: %0.4f\n lam3: %0.4f\n'% (N,ipr1,ipr2,ipr3,lam1,lam2,lam3)                                                                                                       
   print(nx.is_connected(G))
   return G
#eigval_bound(M)
  
