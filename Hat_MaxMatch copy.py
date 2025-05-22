from hopcroftkarp import HopcroftKarp
from networkx.algorithms import bipartite
from scipy.spatial import KDTree as KD
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
import itertools as it
import pickle 
import time
from math import radians, sin, cos


def MaxMatch(G, n, Nodes,x):
	def mapping():
		map = {}
		for j in Nodes:
			flag = True
			while flag == True:
				r = random.randint(0,x)
				if(r not in map.values()):
					map[j]=r
					flag = False
		return map
	
	MM=[]
	for j in range(n):
		map = mapping()
		H = nx.relabel_nodes(G, map)
		
		
		X, Y = bipartite.sets(H)

		BP = {}
		for i in X:
			BP[i]=0
		for i in Y:
			BP[i]=1
		nx.set_node_attributes(H, BP, name="bipartite")
	
		GD={}
		for line in bipartite.generate_edgelist(H, data=False):
			l=line.split(' ')
			if(int(l[0]) not in list(GD.keys())):
				GD[int(l[0])] = {int(l[1])}
			else:
				GD[int(l[0])].add(int(l[1]))


		M = HopcroftKarp(GD).maximum_matching(keys_only=True)
		'''
		M = nx.max_weight_matching(H,maxcardinality = True)
		#M = nx.bipartite.maximum_matching(H)
		#M = nx.bipartite.hopcroft_karp_matching(H)
		'''
		rev_map = dict((v,k) for k,v in map.items())

		ME = []
		for u,v in M.items():
			ME.append((rev_map[u],rev_map[v]))

		pickle.dump(ME, open('H_MaxMatches_3out/MaxMatch_%i.txt'%j,'wb'))    

def rotate(x,y,Op):

	if(len(Op) == 1):
		cx = 0
		cy = 0
	else:
		cx = float(Op[1])
		cy = float(Op[2])

	angle = radians(float(Op[0]))
	
	x_new = cos(angle)*(x - cx) - sin(angle)*(y - cy) + cx
	y_new = sin(angle)*(x - cx) + cos(angle)*(y - cy) + cy
	
	return x_new,y_new

def matrix(x,y,Op):
	x_new = x*float(Op[0])+y*float(Op[2])+float(Op[4])
	y_new = x*float(Op[1])+y*float(Op[3])+float(Op[5])
	
	return x_new,y_new

def translate(x,y,Op):
	x_new = x + float(Op[0])
	y_new = y + float(Op[1])
	
	return x_new,y_new

def Pos_transform(t,pos):
	Op_ = T[t]
	Op = Op_[0]
	if(Op == 'matrix'):
		p=[]
		Op_v = Op_[1]
		for u,v in pos:
			u_,v_ = matrix(u,v,Op_v)
			p.append((round(u_,2), round(v_,2)))
		return p
	elif(Op == 'translate'):
		p=[]
		Op_v = Op_[1]
		for u,v in pos:
			u_,v_ = translate(u,v,Op_v)
			p.append((round(u_,2), round(v_,2)))
		return p
	elif(Op == 'rotate'):
		p=[]
		Op_v = Op_[1]
		for u,v in pos:
			u_,v_ = rotate(u,v,Op_v)
			p.append((round(u_,2), round(v_,2)))
		return p
	else:
		return pos


def Hat(name):
	f = open('%s.svg'%name,'r')
	lines = f.readlines()

	for i in range(len(lines))[::-1]:
		if('defs' in lines[i]):
			a=i
			break

	lines = lines[a+2:]
	transform = {}
	for i in range(len(lines)):
		if('use' in lines[i]):
			k = lines[i].replace(' ','')[7:-2]
			if('matrix' in lines[i+1]):
				k1 = lines[i+1].replace(' ','')[18:-4].split(',')
				transform[k]=['matrix',k1]
			elif('translate' in lines[i+1]):
				k1 = lines[i+1].replace(' ','')[21:-4].split(',')
				transform[k]=['translate',k1]
			elif('rotate' in lines[i+1]):
				k1 = lines[i+1].replace(' ','')[18:-4].split(',')
				transform[k]=['rotate',k1]
			else:
				transform[k]=['']

	polygons = []
	p = []
	for i in lines:
		if('use' in i):
			k = i.replace(' ','')[7:-2]
			p.append(k)
		elif('polygon' in i):
			polygons.append(p.copy())
		elif('</g>' in i):
			p.pop()
	return polygons, transform
	
P,T = Hat('H_output-3')

Pos=[]
for i in P:
	pos = [(-1,-1.7320508),(0,-1.7320508),(1,-1.7320508),(1.5,-0.8660254),(3,-1.7320508),(4.5,-0.8660254),(4,0),(3,0),(3,1.7320508),(1.5,2.5980762),(1,1.7320508),(0,1.7320508),(0,0),(-1.5,-0.8660254)]
	for j in i[::-1]:
		pos = Pos_transform(j,pos)
	Pos.append(pos)

pos={}
edges=[]
n=0
for i in Pos:
	for j in range(len(i)):
		if(bool(pos) == False):
			u,v = i[j]
			u1,v1 = i[j+1]
			pos[n] = (u,v)
			pos[n+1] = (u1,v1)
			if((n,n+1) not in edges and (n+1,n) not in edges):
				edges.append((n,n+1))
			n += 2

		elif(j != len(i)-1):
			u,v = i[j]
			u1,v1 = i[j+1]
			
			t = KD(list(pos.values()))
			a = t.query([(u,v)], distance_upper_bound=0.5)  #k=1
			b = t.query([(u1,v1)], distance_upper_bound=0.5)
			if(a[1][0] >= len(pos.values()) and b[1][0] == len(pos.values())):
				pos[n] = (u,v)
				pos[n+1] = (u1,v1)
				if((n,n+1) not in edges and (n+1,n) not in edges):
					edges.append((n,n+1))
				n += 2
			elif(a[1][0] == len(pos.values()) and b[1][0] < len(pos.values())):
				pos[n] = (u,v)
				if((n,list(pos.keys())[b[1][0]]) not in edges and (list(pos.keys())[b[1][0]],n) not in edges):
					edges.append((n,list(pos.keys())[b[1][0]]))
				n += 1
			elif(a[1][0] < len(pos.values()) and b[1][0] >= len(pos.values())):
				pos[n] = (u1,v1)
				if((list(pos.keys())[a[1][0]],n) not in edges and (n,list(pos.keys())[a[1][0]]) not in edges):
					edges.append((list(pos.keys())[a[1][0]],n))
				n += 1
			else:
				if((list(pos.keys())[a[1][0]],list(pos.keys())[b[1][0]]) not in edges and (list(pos.keys())[b[1][0]],list(pos.keys())[a[1][0]]) not in edges):
					edges.append((list(pos.keys())[a[1][0]],list(pos.keys())[b[1][0]]))
		else:
			u,v = i[j]
			u1,v1 = i[0]
			t = KD(list(pos.values()))
			a = t.query([(u,v)], distance_upper_bound=0.5)      
			b = t.query([(u1,v1)], distance_upper_bound=0.5)  
			if((list(pos.keys())[a[1][0]],list(pos.keys())[b[1][0]]) not in edges and (list(pos.keys())[b[1][0]],list(pos.keys())[a[1][0]]) not in edges): 
				edges.append((list(pos.keys())[a[1][0]],list(pos.keys())[b[1][0]]))
		
G = nx.Graph()
G.add_edges_from(edges)
Nodes = list(G.nodes())
X, Y = bipartite.sets(G)
N = len(Nodes)
print(N,len(X),len(Y))

pickle.dump(G,open('Hout3.txt','wb'))
pickle.dump(pos,open('Pos_Hout3.txt','wb'))
nx.draw(G, pos=pos, node_size=5, with_labels = True, font_size=5)
plt.show()

Efsdefds

MM = MaxMatch(G,20,Nodes,N*2)

c=0
DDL = []
for i in range(10,11):
	M = pickle.load(open('H_MaxMatches_3out/MaxMatch_%i.txt'%i,'rb'))
	l = list(range(10))
	#l.remove(i)
	for j in l:
		r=[]
		E = M.copy()
		M1 = pickle.load(open('H_MaxMatches_3out/MaxMatch_%i.txt'%j,'rb'))
		E1 = M1.copy()
		for u,v in M1:
			if((u,v) in M):
				E.remove((u,v))
				E1.remove((u,v))
				r.append((u,v))
			elif((v,u) in M):
				E.remove((v,u))
				E1.remove((u,v))
				r.append((v,u))

		E.extend(E1)
		DDL.append(E)
		c=c+1
		print(c)
	
		
		EC=[]
		W=[]
		for u,v in list(G.edges()):
			if((u,v) in E or (v,u) in E):
				if((u,v) in M or (v,u) in M):
					EC.append('purple')
					W.append(2)
				elif((u,v) in M1 or (v,u) in M1):
					EC.append('purple')
					W.append(2)
			#elif((u,v) in r or (v,u) in r):
			#	EC.append('red')
			#	W.append(2)
			else:
				EC.append('black')
				W.append(1)
		nx.draw(G, node_size= 0, pos=pos, edge_color=EC, width=W) 
		plt.show()
	
#pickle.dump(DDL, open('DDLs_MaxMatchU1_1M_relabel.txt','wb'))

Sfsdf

Dim_E = []
EF=[]
mp=[]
for i in range(20):
	print(i)
	E=[]
	MN=[]
	
	m = pickle.load(open('H_MaxMatches_3out/MaxMatch_%i.txt'%i,'rb'))
	
	for u,v in m:
		MN.append(u)
		MN.append(v)

	DF  = set(Nodes) - set(MN)  
	print(DF)
	
	Dim_E.append(set(m))
	
	for u,v in G.edges():
		if((u,v) not in m and (v,u) not in m):
			E.append((u,v))
	EF.append(set(E))

	'''
	NC = []
	for I in G.nodes():
		if(I in DF):
			NC.append('red')
		else:
			NC.append('black')
	EC = []
	W1=[]
	for u,v in G.edges():
		if((u,v) in m):
			EC.append('purple')
			W1.append(2)
		elif((v,u) in m):
			EC.append('purple')
			W1.append(2)
		else:
			EC.append('black')
			W1.append(1)
			
	nx.draw(G, pos=pos, node_size=10, edge_color=EC, width=W1, with_labels=True,font_size=5, node_color = NC)
	plt.show()
	'''

com_EF = set.intersection(*Dim_E)
com_EF1 = set.intersection(*EF)
com = com_EF.copy()
com.update(com_EF1)
Rem = []
for u,v in G.edges():
	if((u,v) not in com and (v,u) not in com):
			Rem.append((u,v))

EC = []
W1 = []
for u,v in G.edges():
	if((u,v) in com_EF):
		EC.append('purple')
		W1.append(2)
	elif((v,u) in com_EF):
		EC.append('purple')
		W1.append(2)
	elif((u,v) in com_EF1):
		EC.append('green')
		W1.append(2)
	elif((v,u) in com_EF1):
		EC.append('green')
		W1.append(2)
	elif((u,v) in Rem):
		EC.append('red')
		W1.append(2)
	elif((v,u) in Rem):
		EC.append('red')
		W1.append(2)
	else:
		EC.append('black')
		W1.append(0.7)
	
			
nx.draw(G, pos=pos, edge_color=EC, width=W1, node_size=0)
plt.show()





