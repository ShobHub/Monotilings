from hopcroftkarp import HopcroftKarp
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
import itertools as it
import pickle 
import time

import sys
sys.setrecursionlimit(3000)

'''
def Spectre(name):
	Pos = {}
	n=0
	edges=[]

	f = open("%s.svg"%name, "r")
	lines = f.readlines()
	c=0
	for i in lines[5:]:
		l = i.split()
		if(l[0] == '<polygon'): # and len(l[1:-3])==13):
			#print(c)
			c=c+1
			l = l[1:-3]
			l[0] = l[0][8:]
			#print(len(l))
			
			for j in range(len(l[:-1])):
				x,y = l[j].split(',')
				x1,y1 = l[j+1].split(',')
				if('"' in y1):
					y1 = y1[:-1]
				u,v = round(float(x),2), round(float(y),2)
				u1,v1 = round(float(x1),2), round(float(y1),2)

				if((u,v) not in Pos.values()):
					Pos[n] = (u,v)
					a=n
					n=n+1
				
				else:
					a = list(Pos.keys())[list(Pos.values()).index( (u,v) )]

				if((u1,v1) not in Pos.values()):
					Pos[n] = (u1,v1)
					b=n
					n=n+1
				else:
					b = list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
				edges.append((a,b))

			x,y = l[-1].split(',')
			x1,y1 = l[0].split(',')
			u,v = round(float(x),2), round(float(y[:-1]),2)
			u1,v1 = round(float(x1),2), round(float(y1),2)
			a = list(Pos.keys())[list(Pos.values()).index( (u,v) )]
			b = list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
			edges.append((a,b))
			

	#print(Pos)
	
	Pos1=Pos.copy()
	for u in list(Pos1.keys())[13:]: 
		x,y = Pos1[u]
		x_new = -0.24999999999999997 * x + 0.4330127018922194 * y + 2.5
		y_new = 0.4330127018922194 * x + 0.24999999999999997 * y + 0.8660254037844386
		
		Pos[u] = (x_new, y_new)
	
	E=[]
	E = edges.copy()
	for u,v in edges:
		E.append((u+13,v+13))
		x_new = 1*Pos[u][0]+0*Pos[u][1]+3
		y_new = 0*Pos[u][0]-1*Pos[u][1]+1.7255
		x_new1 = 1*Pos[v][0]+0*Pos[v][1]+3
		y_new1 = 0*Pos[v][0]-1*Pos[v][1]+1.7255

		Pos[u+13] = (x_new, y_new)
		Pos[v+13] = (x_new1, y_new1)

	G = nx.Graph()
	G.add_edges_from(edges)
	
	return G, Pos
'''

def Spectre(name):

	f = open("%s.svg"%name, "r")
	lines = f.readlines()
	D={}
	c=0

	while c<len(lines):
		l = lines[c].split()
		if(l[0]=='<polygon'):
			a=c
			D[a]=[]
		elif(l[0]=='<use'):
			ex=[]
			ex.append(float(l[2][18:]))
			ex.append(float(l[3]))
			ex.append(float(l[4]))
			ex.append(float(l[5]))
			ex.append(float(l[6]))
			ex.append(float(l[7][:-4]))
			D[a].append(ex)
		c=c+1
	
	D[list(D.keys())[-1]].pop()
	return D, lines

def Transform(x,y,Op):
	x_new = x*Op[0]+y*Op[2]+Op[4]
	y_new = x*Op[1]+y*Op[3]+Op[5]
	
	return x_new,y_new
			 

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
		
		#M = nx.bipartite.maximum_matching(H)
		#M = nx.bipartite.hopcroft_karp_matching(H)
		
		rev_map = dict((v,k) for k,v in map.items())

		ME = []
		for u,v in M.items():
			ME.append((rev_map[u],rev_map[v]))
	
		#MM.append(ME)	
		pickle.dump(ME, open('Hat_MaxMatches_2out/MaxMatch_%i.txt'%j,'wb'))

	#return MM

'''
DFL=[]
NL=[]
for i in range(1,3):
	
	print(i)
	
	G, Pos = Spectre('Hat_output-%i'%i)
	Nodes = list(G.nodes())
	N = len(Nodes)
	print(N)
	NL.append(N)
	
	MM = MaxMatch(G,10,Nodes,N*2)

	s=0
	for j in MM:
		MN=[]
		for u,v in j:
			MN.append(u)
			MN.append(v)

		DF  = set(Nodes) - set(MN)    
		s=s+len(DF)
	
	DFL.append((s/10)/N)
	

pickle.dump([NL,DFL], open('Hat_MonomerDensity_avg10.txt','wb'))
plt.scatter(NL, DFL, c='steelblue')
plt.show()
'''

D, lines = Spectre('H_output-2')
Polygons = []
PolyVisited = []
for i in D.keys():
	l = lines[i].split()
	l = l[1:-3]
	l[0] = l[0][8:]
	if(l not in PolyVisited):
		PolyVisited.append(l)
		p=[]
		for j in range(len(l)):
			x,y = l[j].split(',')
			if('"' in y):
				y = y[:-1]
			u,v = round(float(x),2), round(float(y),2)
			p.append((u,v))
		print(len(p))
		Polygons.append(p)

for i in D.keys():
	if(D[i]!=[]):
		l = lines[i].split()
		l = l[1:-3]
		l[0] = l[0][8:]
		p=[]
		for j in range(len(l)):
			x,y = l[j].split(',')
			if('"' in y):
				y = y[:-1]
			u,v = round(float(x),2), round(float(y),2)
			p.append((u,v))	
		for k in D[i]:
			p1=[]
			for u,v in p:
				un,vn = Transform(u,v,k)
				p1.append((round(float(un),2), round(float(vn),2)))
			Polygons.append(p1)

#print(Polygons)

Pos={}
n=0
edges=[]

kkk=[[1,0, 0, 1, -4.499999999999998, -2.5980762113533156],[1, 2.7755575615628914e-16, -2.7755575615628914e-16, 1, 1.7763568394002505e-15, -3.4641016151377553],[-0.49999999999999994, -0.8660254037844388, 0.8660254037844388, -0.49999999999999994, -3.552713678800501e-15, 5.196152422706632]]

kkk1=[[-0.5, -0.8660254037844388, 0.8660254037844388, -0.5, -2.0000000000000036, 5.196152422706632],[-0.5, -0.8660254037844386, 0.8660254037844386, -0.5, -4.999999999999998, 1.7320508075688772],[1, 2.9605947323337506e-16, -2.9605947323337506e-16, 1, 1.0000000000000018, -5.196152422706632],[-0.5000000000000001, 0.8660254037844387, -0.8660254037844387, -0.5000000000000001, 3.9999999999999964, 3.464101615137754],[0.5000000000000004, -0.8660254037844384, 0.8660254037844384, 0.5000000000000004, 2.5, 0.8660254037844379],[1, 0, 0, 1, -3.5, -4.332]]

for i in kkk:
	for j in [9,10,11,12]:
		p1=[]
		for u,v in Polygons[j]:
			un,vn = Transform(u,v,i)
			p1.append((round(float(un),2), round(float(vn),2)))
		Polygons.append(p1)

for i in kkk1:
	for j in [19,20]:
		p1=[]
		for u,v in Polygons[j]:
			un,vn = Transform(u,v,i)
			p1.append((round(float(un),2), round(float(vn),2)))
		Polygons.append(p1)


ran = it.chain(range(19,19),range(-24,0))
for i in ran:   #range(len(Polygons)): 
	j = Polygons[i]
	for k in range(len(j)-1):
		u,v = j[k]
		u1,v1 = j[k+1]

		if((u,v) not in Pos.values()):
			Pos[n] = (u,v)
			a=n
			n=n+1
				
		else:
			a = list(Pos.keys())[list(Pos.values()).index( (u,v) )]

		if((u1,v1) not in Pos.values()):
			Pos[n] = (u1,v1)
			b=n
			n=n+1
		else:
			b = list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
		edges.append((a,b))

	x,y = j[-1]
	x1,y1 = j[0]
	a = list(Pos.keys())[list(Pos.values()).index( (x,y) )]
	b = list(Pos.keys())[list(Pos.values()).index( (x1,y1) )]
	edges.append((a,b))

G=nx.Graph()
G.add_edges_from(edges)
print(len(G.nodes()))
nx.draw(G,pos=Pos,node_size=0) #,with_labels=True)
plt.show()

dxfdds


#MM = MaxMatch(G,100,Nodes,N*2)

'''
c=0
DDL = []
for i in range(10,11):
	M = pickle.load(open('Hat_MaxMatches_2out/MaxMatch_%i.txt'%i,'rb'))
	l = list(range(90,100))
	#l.remove(i)
	for j in l:
		E = M.copy()
		M1 = pickle.load(open('Hat_MaxMatches_2out/MaxMatch_%i.txt'%j,'rb'))
		E1 = M1.copy()
		for u,v in M1:
			if((u,v) in M):
				E.remove((u,v))
				E1.remove((u,v))
			elif((v,u) in M):
				E.remove((v,u))
				E1.remove((u,v))

		E.extend(E1)
		DDL.append(E)
		c=c+1
		print(c)
	
		
		EC=[]
		W=[]
		for u,v in list(G.edges()):
			if((u,v) in E or (v,u) in E):
				EC.append('purple')
				W.append(2)
			else:
				EC.append('black')
				W.append(1)
		nx.draw(G, node_size= 0, pos=Pos, edge_color=EC, width=W) 
		plt.show()
	
#pickle.dump(DDL, open('DDLs_MaxMatchU1_1M_relabel.txt','wb'))
'''

Dim_E = []
EF=[]
for i in range(100):
	print(i)
	E=[]
	MN=[]
	
	m = pickle.load(open('Hat_MaxMatches_2out/MaxMatch_%i.txt'%i,'rb'))
	
	Dim_E.append(set(m))

	for u,v in G.edges():
		if((u,v) not in m and (v,u) not in m):
			E.append((u,v))
	EF.append(set(E))

	'''
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
			
	nx.draw(G, pos=Pos, node_size=0, edge_color=EC, width=W1)
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
	
			
nx.draw(G, pos=Pos, edge_color=EC, width=W1, node_size=0)
plt.show()

'''
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

'''
































#############  Alternate ####################
'''
st = time.time()
	for i in lines[2:-2:2]:
		l = i.split()[1:-3]
		l[0] = l[0][8:]
		for j in l:
			x,y = j.split(',')
			if('"' in y):
				y=y[:-1]
			u,v = round(float(x),2), round(float(y),2)
			if((u,v) not in Pos.values()):
				Pos[n] = (u,v)
				n=n+1

	for i in lines[2:-2:2]:
		l = i.split()[1:-3]
		l[0] = l[0][8:]
		for j in range(len(l[:-1])):
			x,y = l[j].split(',')
			x1,y1 = l[j+1].split(',')
			if('"' in y1):
				y1 = y1[:-1]
			u,v = round(float(x),2), round(float(y),2)
			u1,v1 = round(float(x1),2), round(float(y1),2)
			a = list(Pos.keys())[list(Pos.values()).index( (u,v) )]
			b = list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
			edges.append((a,b))
		x,y = l[-1].split(',')
		x1,y1 = l[0].split(',')
		u,v = round(float(x),2), round(float(y[:-1]),2)
		u1,v1 = round(float(x1),2), round(float(y1),2)
		a = list(Pos.keys())[list(Pos.values()).index( (u,v) )]
		b = list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
		edges.append((a,b))

	print(time.time()-st)
'''
	
