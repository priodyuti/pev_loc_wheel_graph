# email for communication: prioyutipradhan@gmail.com
import networkx as nx                                              
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import codecs
import sys
from operator import itemgetter
from numpy import linalg as LA
from functions_base import store_undirected_network, read_network

def get_ipr(A):
  e_val, e_vec = LA.eigh(A)
  d = np.shape(e_vec)   
  n = d[0]   
  lam1 = e_val[n-1]  
  Ev1 = e_vec[:,n-1]
  
          
  IPR1 = 0.0;
  for i in range(n):
    IPR1 = IPR1 + pow(Ev1[i],4);
  del e_val
  del e_vec
  del Ev1
 
  return IPR1,lam1


path = 'wheel_random_comb_loc.txt'

s1 = 398
r = 11
Flag = False
while not Flag:
  G1 = nx.random_regular_graph(r,s1)
  #G1 = nx.fast_gnp_random_graph(s1, p)
  Flag = nx.is_connected(G1)
  print Flag

M1 = np.zeros((s1,s1),dtype=int)
Edge_list1 = nx.edges(G1)
for edge in Edge_list1:
  x = edge[0]
  y = edge[1]
  M1[x-1][y-1] = 1                                                                        
  M1[y-1][x-1] = 1

ipr11,lam11 = get_ipr(M1)

s2 = 100

G2 = nx.wheel_graph(s2)
   
N = s1 + s2 + 1
   
M = np.zeros((N,N), dtype=int)
M2 = np.zeros((s2,s2), dtype=int)

Edge_list2 = nx.edges(G2)
for edge in Edge_list2:
     x = edge[0]
     y = edge[1]
     M2[x][y] = 1
     M2[y][x] = 1
  
   #print 'Number of nodes in wheel graph: %d' %s2
   #print 'Number of edges in wheel graph: ', nx.number_of_edges(G)
   
ipr21,lam21 = get_ipr(M2)
   #print 's2: %d lam11: %0.4f lam12: %0.4f lam21: %0.4f lam22: %0.4f\n'%(s2,lam11,lam12,lam21,lam22)                
   #print 's2: %d ipr1: %0.4f ipr2: %0.4f ipr3: %0.4f lam1: %0.4f lam2: %0.4f lam3: %0.4f\n'% (s2,ipr1,ipr2,ipr3,lam1,lam2,lam3)                                                                                                        

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

gp = nx.from_numpy_matrix(M)
ipr1,lam1 = get_ipr(M)
print(nx.is_connected(gp))
   
print 'N: %d s2: %d ipr1: %0.4f lam1: %0.4f lam11: %0.4f lam21: %0.4f\n'%  (N,s2,ipr1,lam1,lam11,lam21)                
#fd.write(str(N)+' '+str(s2)+' '+str(ipr1)+' '+str(lam1)+' '+str(lam11)+' '+str(lam21)+'\n')            
store_undirected_network(M, N, path)

del G2
del gp
del M
del M2   
   
#fd.close()


