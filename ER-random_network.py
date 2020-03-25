# email for communication: prioyutipradhan@gmail.com
import networkx as nx                                              
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import codecs
import sys
from operator import itemgetter
from numpy import linalg as LA
from functions_base import store_undirected_network


n = 1524
k = 6
p = k/float(n)

Flag = False
while not Flag:
  G = nx.fast_gnp_random_graph(n, p)
  Flag = nx.is_connected(G)
  print Flag

Edges = nx.number_of_edges(G)
avg_deg = Edges*2/float(n)
print 'Number of nodes in ER random graph: %d' %n
print 'Number of edges in ER random graph: ', nx.number_of_edges(G)
print 'Average degree of ER random graph: %f' %avg_deg

M = np.zeros((n,n), dtype=int)

Edge_list = nx.edges(G)
for edge in Edge_list:
  x = edge[0]
  y = edge[1]
  M[x][y] = 1
  M[y][x] = 1

store_undirected_network(M, n, 'ER_random.txt')

