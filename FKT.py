import matplotlib.pyplot as plt
from termcolor import colored
from scipy.spatial import KDTree
import scipy as spp
import sympy as sym
import networkx as nx
import numpy as np
import pickle 
import random
import time

def FKT(G,Pos,d_pos,PL,BDL): #n):
	'''
 	G: Lattice Graph
  	Pos: Nodes positions dictionary
   	d_pos: Dual lattice nodes position dictionary
    	PL: List of Plaquettes
     
 	Implement High-level description of the FKT algorithm steps mentioned here: https://en.wikipedia.org/wiki/FKT_algorithm
  	'''
	St = time.time()
	# G = pickle.load(open('Delta/output-%i.pickle'%n,'rb'))
	# Pos = pickle.load(open('Delta/Pos_output-%i.pickle'%n,'rb'))
	# d_pos = pickle.load(open('Delta/dual_pos_output-%i.pickle'%n,'rb'))
	# PL = pickle.load(open('Delta/Plaqs_output-%i.pickle'%n,'rb'))
	# BDL = pickle.load(open('Delta/BoundaryLoop_output-%i.txt'%n,'rb'))
	
	#G = nx.DiGraph()
	s=0
	e=[]
	for i in range(len(BDL)-1):
		s = s + np.array(Pos[BDL[i]])
		e.append((BDL[i],BDL[i+1]))
	e.append((BDL[-1],BDL[0]))
	s = (s + np.array(Pos[BDL[-1]]))/len(BDL)
	PL['out'] = e
	d_pos['out'] = s  #(1069,1069)     
	
	T1 = nx.minimum_spanning_tree(G)
	T1E = list(T1.edges())
	
	T2 = nx.Graph()
	T2.add_nodes_from(d_pos.keys())
	kdtree = KDTree(list(d_pos.values()))                                        
	dk = list(d_pos.keys())

	def com_edges(a,b):
		CE=[]
		for u,v in a:
			if((u,v) in b or (v,u) in b):
				CE.append((u,v))
		return CE

	def NinT1(CE):
		flag = 0
		for u,v in CE:
			if((u,v) in T1.edges() or (v,u) in T1.edges()):
				flag = 0
			else:
				flag = 1
				break
		if(flag == 1):
			return True
		else:
			return False

	for i in T2.nodes():
		sp = PL[i]
		xu,xv = d_pos[i]
		dis, ind = kdtree.query([xu, xv], k=240, distance_upper_bound=100)   #36  240(RG 884)
		ind = list(filter(lambda num: num < len(T2.nodes()), ind)) 
		#T2N = list(T2.nodes())
		#T2N.remove(i)
		for j in ind[1:]:     #T2N:
			ep = PL[dk[j]]        #j
			CE = com_edges(sp,ep)
			if(CE!=[] and NinT1(CE) == True):
				T2.add_edge(i,dk[j])         #j
	
	flag = False
	for I in T2.nodes():
		ep = PL[I]   
		for u,v in ep:
			if(((u,v) in e or (v,u) in e) and  NinT1([(u,v)]) == True and I!='out'):
				flag = True
				T2.add_edge(I,'out')
				break
		# if flag:
		# 	break

	# nx.draw(G,pos=Pos,node_size=0,width=0.7)
	# nx.draw(T1,pos=Pos,node_size=0,edge_color='red',width=1)
	# nx.draw(T2,pos=d_pos,node_size=5, node_color='blue', edge_color='blue', width=1)  #, with_labels=True, font_size=8)
	# plt.show()

	def Orien(a,b,p):
		vec1 = np.array(Pos[a]) - np.array(p)
		vec2 = np.array(Pos[b]) - np.array(p)
		CP = vec1[0]*vec2[1]-vec1[1]*vec2[0]
		if(CP>0):
			return 'AC'
		if(CP<0):
			return 'C'

	# T2N = list(T2.nodes())
	# T2N.remove('out')
	print('&')
	while (list(T2.nodes()) != ['out']):  #1): 
		# print(len(T2.nodes()))  #T2.nodes(),
		Nd1 = list(filter(lambda num: (T2.degree(num) == 1 and num!='out'), T2.nodes())) 
		#print(len(Nd1))
		#Nd1 = []
		#for i in T2.nodes():
		#	if(T2.degree(i) == 1 and i!='out'):
		#		Nd1.append(i)
		#cc = 0
		for i in Nd1:
			#print(cc)
			#cc += 1
			s_C = 0
			s_AC = 0
			nT1 = {}
			for u,v in PL[i]:
				if((u,v) in T1E):
					d = Orien(u,v, d_pos[i])
					if(d == 'C'):
						s_C = s_C+1
					else:
						s_AC = s_AC+1
				elif((v,u) in T1E):
					d = Orien(v,u, d_pos[i])
					if(d == 'C'):
						s_C = s_C+1
					else:
						s_AC = s_AC+1
				else:
					nT1[(u,v)] = Orien(u,v, d_pos[i])
			#print(i,s_C,s_AC)
			if(s_C%2 != 0):
				for u,v in nT1.keys():
					if(nT1[(u,v)] == 'C'):
						T1E.append((v,u))
					else:
						T1E.append((u,v))
			else:
				flag = True
				if('C' in nT1.values()):
					for u,v in nT1.keys():
						if(nT1[(u,v)] == 'C' and flag == True):
							T1E.append((u,v))
							flag = False
						elif(nT1[(u,v)] == 'C' and flag == False):
							T1E.append((v,u))
						elif(nT1[(u,v)] == 'AC'):
							T1E.append((u,v))
				else:
					flag = True
					for u,v in nT1.keys():
						if(nT1[(u,v)] == 'AC' and flag == True):
							T1E.append((v,u))
							flag = False
						else:
							T1E.append((u,v))
				
			T2.remove_node(i)
	'''
	##################### TEST ##############################################
	for i in T2N:
		s_C=0
		s_AC=0
		for u,v in PL[i]:
			if((u,v) in T1E):
				d = Orien(u,v, d_pos[i])
				if(d == 'C'):
					s_C = s_C+1
					print(u,v,'C')
				else:
					s_AC = s_AC+1
					print(u,v,'AC')
			elif((v,u) in T1E):
				d = Orien(v,u, d_pos[i])
				if(d == 'C'):
					s_C = s_C+1
					print(v,u, 'C')
				else:
					s_AC = s_AC+1
					print(v,u, 'AC')
		print(i,s_C,s_AC)
		if(s_C%2 == 0):
			print(colored('*','red'))
	##################### TEST ##############################################
	'''
	print('&')
	g = nx.DiGraph()
	g.add_edges_from(T1E)
	#pickle.dump(g, open('FKT_graph_6.txt','wb'))
	#M_set = set(map(tuple, sorted(T1E)))
	#for u,v in G.edges():
	#	if ((u,v) not in M_set and (v,u) not in M_set):
	#		print(u,v)
	# nx.draw(g,pos=Pos,node_size=2, with_labels=True,font_size=8)
	# plt.show()

	A = nx.to_numpy_array(g)
	A = A - A.T
	det = np.linalg.det(A)
	pfa = abs(np.sqrt(det))
	# f = sym.log(pfa) * (2/len(g.nodes()))
	# W = sym.exp(f)

	return pfa  #f,W,len(g.nodes())

# fe = []
# Wm = []
# Nn = []
'''
for i in range(1,5):
	f,W,N  = FKT(i)
	print(f,W,N)
	#fe.append(f)
	#Wm.append(W)
	#Nn.append(N)
'''
# g = pickle.load(open('FKT_graph_6.txt','rb'))
# A = nx.to_scipy_sparse_array(g,format='csr')
# A = A-A.transpose()
# Lu = spp.sparse.linalg.splu(A)
# Det = sym.prod(sym.Matrix(np.diag(Lu.U.A)))
# pfa = abs(sym.sqrt(Det))
# print(pfa)

# f = sym.log(pfa) * (2/len(g.nodes()))
# print(f)
# Sdfs

# W = sym.exp(f)
# fe.append(f)
# Wm.append(W)
# Nn.append(214662)
# print(f,W, 214662)

#pickle.dump([Nn,fe],open('FreeEnergy.txt','wb'))
#pickle.dump([Nn,Wm],open('MolecularFreedom.txt','wb'))

#Nn,fe = pickle.load(open('FreeEnergy.txt','rb'))
#Nn,Wm = pickle.load(open('MolecularFreedom.txt','rb'))

# plt.scatter(Nn,fe)
# plt.plot(Nn,fe)
# plt.xlabel('Number of Nodes')
# plt.ylabel('Free energy per dimer (f)')
# #plt.ylim(0.02,0.04)
# plt.show()
# plt.scatter(Nn,Wm)
# plt.plot(Nn,Wm)
# plt.xlabel('Number of Nodes')
# plt.ylabel('Molecular freedom (W)')
#plt.ylim(1.015,1.05)
# plt.show()


'''
g = nx.DiGraph()
g.add_edges_from([(1,9),(9,8),(3,8),(11,3),(2,11),(2,1),(3,10),(10,5),(5,8),(4,2),(4,6),(6,12),(12,5),(4,13),(13,7),(7,6),(7,14),(1,14)])
print(len(g.edges()))
A = np.zeros((14,14))
for u,v in g.edges():
	A[u-1,v-1] = 1
print(A)
Det = np.linalg.det(A-A.T)
Pfa = abs(np.sqrt(Det))
print(Pfa)
'''

	






