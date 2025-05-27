from hopcroftkarp import HopcroftKarp
import xml.etree.ElementTree as ET
from scipy.spatial import cKDTree
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
from math import radians, sin, cos
from scipy.spatial import KDTree
import networkx as nx
import numpy as np
import random
import pickle 
import time
import pylab
from matplotlib.pyplot import pause

import sys
sys.setrecursionlimit(10000)

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

def Pos_transform(t,u,v):
	Op = t[0]
	if(Op == 'matrix'):
		Op_v = t[1]
		u_,v_ = matrix(u,v,Op_v)
		return(int(u_),int(v_))       #round(u_,1), round(v_,1))
	elif(Op == 'translate'):
		Op_v = t[1]
		u_,v_ = translate(u,v,Op_v)
		return(int(u_),int(v_))       #(round(u_,1), round(v_,1))
	elif(Op == 'rotate'):
		Op_v = t[1]
		u_,v_ = rotate(u,v,Op_v)
		return(int(u_),int(v_))       #(round(u_,1), round(v_,1))
	else:
		return u,v
'''
def Spectre(name):
	Pos = {}
	n=0
	edges=[]
	c=0
	f = open("%s.svg"%name, "r")
	lines = f.readlines()
	LL = len(lines)    #*
	d_pos = {}
	PL = {}
	for i in range(35,LL-2,7):    #i in lines[2:-2:2]:   
		print(i)
		l = lines[i].split()[:-1]    #i.split()[1:-3]
		l[0] = l[0][8:]
		
		if('matrix' in lines[i+5]):                                         #*******
			k11 = lines[i+5].replace(' ','')[18:-5].split(',')
			id = ['matrix',k11]
		elif('translate' in lines[i+5]):
			k11 = lines[i+5].replace(' ','')[21:-5].split(',')
			id = ['translate',k11]
		elif('rotate' in lines[i+5]):
			k11 = lines[i+5].replace(' ','')[18:-5].split(',')
			id = ['rotate',k11]                                         #*******
		
		EP = []
		s = np.array((0,0))
		for j in range(len(l[:-1])):
			x,y = l[j].split(',')
			x1,y1 = l[j+1].split(',')
			if('"' in y1):
				y1 = y1[:-1]
			#u,v = round(float(x),2), round(float(y),2)
			#u1,v1 = round(float(x1),2), round(float(y1),2)
			u,v = Pos_transform(id, round(float(x),1), round(float(y),1))
			u1,v1 = Pos_transform(id, round(float(x1),1), round(float(y1),1))
			
			s = s + np.array((u,v))
			if len(Pos.values())!=0:
				kdtree = KDTree(list(Pos.values()))                                         #*******
				dis, ind = kdtree.query([u, v], k=1, distance_upper_bound=3)                
			else:
				ind = 1                                                       #*******
			if(ind >= len(Pos.values())):                                         #(u,v) not in Pos.values()):
				Pos[n] = (u,v)
				a=n
				n=n+1
				
			else:
				a = list(Pos.keys())[ind]          #list(Pos.keys())[list(Pos.values()).index( (u,v) )]

			kdtree = KDTree(list(Pos.values()))                                   #*******
			dis, ind = kdtree.query([u1, v1], k=1, distance_upper_bound=3)      #*******
			if(ind >= len(Pos.values())):                                       #(u1,v1) not in Pos.values()):
				Pos[n] = (u1,v1)
				b=n
				n=n+1
			else:
				b = list(Pos.keys())[ind]        #list(Pos.keys())[list(Pos.values()).index( (u1,v1) )]
			edges.append((a,b))
			EP.append((a,b))

		x,y = l[-1].split(',')
		x1,y1 = l[0].split(',')
		#u,v = round(float(x),2), round(float(y[:-1]),2)     
		#u1,v1 = round(float(x1),2), round(float(y1),2)
		u,v = Pos_transform(id, round(float(x),1), round(float(y),1))
		u1,v1 = Pos_transform(id, round(float(x1),1), round(float(y1),1))
		dis, ind = kdtree.query([u, v], k=1, distance_upper_bound=3)        #*******
		a = list(Pos.keys())[ind]                                          #[list(Pos.values()).index( (u,v) )]
		dis, ind = kdtree.query([u1, v1], k=1, distance_upper_bound=3)      #*******
		b = list(Pos.keys())[ind]                                        #[list(Pos.values()).index( (u1,v1) )]
		edges.append((a,b))
		EP.append((a,b))
		
		s = s + np.array((u,v))
		d_pos[c] = s/14 
		PL[c] = EP
		c=c+1

	G = nx.Graph()
	G.add_edges_from(edges)
	
	return G, Pos, d_pos, PL
'''

def MaxMatch(G,n,Nodes,x):
	''' Generates different maximum matchings of any lattice by relabeling nodes'''
	
	#NO = Nodes.copy()
	def mapping():
		#random.shuffle(NO)
		#map1 = dict(map(lambda i,j : (i,j), Nodes, NO))
		map = {}
		# for j in Nodes:
		# 	flag = True
		# 	while flag == True:
		# 		r = random.randint(0,x)
		# 		if(r not in map.values()):
		# 			map[j]=r
		# 			flag = False
		for i in G.edges():  #Edges:
			map[i] = random.uniform(0.1,1.1)
		nx.set_edge_attributes(G, values=map, name='weight')
		return map
	MM=[]
	for j in range(n):
		mapping()
		# mapp = mapping()
		# H = nx.relabel_nodes(G, mapp)
		# X, Y = bipartite.sets(H)
		# print(len(X),len(Y))
		# BP = dict(map(lambda i : (i,0), X)) | dict(map(lambda i : (i,1), Y))
		# BP = {}
		# for i in X:
		# 	BP[i]=0
		# for i in Y:
		# 	BP[i]=1
		# nx.set_node_attributes(H, BP, name="bipartite")
		# GD={}
		# for line in bipartite.generate_edgelist(H, data=False):
		# 	l=line.split(' ')
		# 	if(int(l[0]) not in list(GD.keys())):
		# 		GD[int(l[0])] = {int(l[1])}
		# 	else:
		# 		GD[int(l[0])].add(int(l[1]))
		# M = HopcroftKarp(GD).maximum_matching(keys_only=True)
		# M = nx.bipartite.maximum_matching(H)
		# M = nx.bipartite.hopcroft_karp_matching(H)
		M = nx.max_weight_matching(G,maxcardinality=True,weight='weight')
		# rev_map = dict((v,k) for k,v in mapp.items())
		ME = []
		for u,v in M:  #.items():
			ME.append(tuple(sorted((u,v))))    #(rev_map[u],rev_map[v]))
		# MM.append(ME)
		pickle.dump(ME, open('Delta/MaxMatches_3out_NB/MaxMatch_%i.txt'%j,'wb'))
	return MM
'''
DFL=[]
NL=[]
for i in range(1,5):
	print(i)
	G, Pos = Spectre('SpectreGamma_output-%i'%i)
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
pickle.dump([NL,DFL], open('SpectreGamma_MonomerDensity_avg10.txt','wb'))
plt.scatter(NL, DFL, c='steelblue')
plt.show()
'''
def parse_svg(svg_file):
    polygons = []
    tree = ET.parse(svg_file)
    root = tree.getroot()
    for polygon in root.findall(".//{http://www.w3.org/2000/svg}polygon"):
        points = polygon.get("points").split()
        polygon_points = [tuple(map(float, point.split(","))) for point in points]
        polygons.append(polygon_points)
    return polygons

def Spectre1(f):
	'''Create Spectre Lattice from the SVG data files obtained from https://cs.uwaterloo.ca/~csk/spectre/ 
 
 	Returns 
  	G: networkx graph
  	Pos: Dictionary for nodes labels and positions
   	d_pos: Dictionary for dual lattice nodes labels and positions
    	PL: List of plaquettes 
	'''
	
	Pos = {}
	n=0
	c=0
	edges=[]
	d_pos = {}
	PL = {}
	polygons = parse_svg("%s.svg"%f)
	print(len(polygons))
	for l in polygons:
		print(c)
		EP = []
		s = np.array((0,0))
		for j in range(len(l[:-1])):
			x,y = l[j]
			x1,y1 = l[j+1]
			u,v = round(float(x),2), round(float(y),2)
			u1,v1 = round(float(x1),2), round(float(y1),2)

			s = s + np.array((u,v))
			if len(Pos.values())!=0:
				kdtree = KDTree(list(Pos.values()))
				dis, ind = kdtree.query([u, v], k=1, distance_upper_bound=3)
			else:
				ind = 1
			if(ind >= len(Pos.values())):
				Pos[n] = (u,v)
				a=n
				n=n+1
			else:
				a = list(Pos.keys())[ind]
			kdtree = KDTree(list(Pos.values()))
			dis, ind = kdtree.query([u1, v1], k=1, distance_upper_bound=3)
			if(ind >= len(Pos.values())):
				Pos[n] = (u1,v1)
				b=n
				n=n+1
			else:
				b = list(Pos.keys())[ind]
			edges.append((a,b))
			EP.append((a,b))

		x,y = l[-1]
		x1,y1 = l[0]
		u,v = round(float(x),2), round(float(y),2)
		u1,v1 = round(float(x1),2), round(float(y1),2)
		dis, ind = kdtree.query([u, v], k=1, distance_upper_bound=3)
		a = list(Pos.keys())[ind]
		dis, ind = kdtree.query([u1, v1], k=1, distance_upper_bound=3)
		b = list(Pos.keys())[ind]
		edges.append((a,b))
		EP.append((a,b))

		s = s + np.array((u,v))
		d_pos[c] = s/14
		PL[c] = EP
		c=c+1

	G = nx.Graph()
	G.add_edges_from(edges)

	return G, Pos, d_pos, PL
'''
G, Pos, d_pos, PL = Spectre1('Delta/output-3')
pickle.dump(G, open('Delta/output-600.pickle', 'wb'))
pickle.dump(Pos, open('Delta/Pos_output-600.pickle', 'wb'))
pickle.dump(d_pos, open('Delta/dual_pos_output-600.pickle', 'wb'))
pickle.dump(PL, open('Delta/Plaqs_output-600.pickle', 'wb'))
'''
G = pickle.load(open('Delta/output-3.pickle', 'rb'))
Pos = pickle.load(open('Delta/Pos_output-3.pickle', 'rb'))
BL = pickle.load(open('Delta/BoundaryLoop_output-3.txt', 'rb'))
#d_pos = pickle.load(open('Delta/dual_pos_output-3.pickle', 'rb'))
#PL = pickle.load(open('Delta/Plaqs_output-3.pickle', 'rb'))

Nodes = list(G.nodes())
Edges = list(G.edges())
N = len(Nodes)
print(N)
n = []
e = []
# for i in range(1,len(BL)-1):
# 	a = np.array(Pos[BL[i-1]]) - np.array(Pos[BL[i]])
# 	b = np.array(Pos[BL[i+1]]) - np.array(Pos[BL[i]])
# 	if(round(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)),1) == -1.0 and BL[i]!=908):
# 		n.append(BL[i])
# 		e.append((BL[i-1],BL[i+1]))
for i in Nodes:
	if(G.degree(i)==2 and i!=908):
		nn = list(G.neighbors(i))
		a = np.array(Pos[nn[0]]) - np.array(Pos[i])
		b = np.array(Pos[nn[1]]) - np.array(Pos[i])
		if (round(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), 1) == -1.0):
			n.append(i)
			e.append((nn[0],nn[1]))
G.remove_nodes_from(n)
G.add_edges_from(e)
Nodes = list(G.nodes())
N = len(Nodes)
print(N,nx.is_bipartite(G))

'''
G = pickle.load(open('Delta/output-6.pickle', 'rb'))
res = dict(map(lambda i,j : (i,j) , list(G.nodes()),list(range(1665406)) ))
G = nx.relabel_nodes(G, res)
pickle.dump(G, open('Delta/output-6_fkt.pickle', 'wb'))

res = dict(map(lambda i,j : ((int(i[0]),int(i[1])),j) , list(G.nodes()),list(range(1665406)) ))
Pos = dict((v,u) for u,v in res.items())

import svgwrite
def create_dual_graph(G, polygons):
    dual_graph = nx.Graph()
    dual_node_positions = {}  # Dictionary to store dual node labels and their positions
    dual_node_edges = {}      # Dictionary to store dual node labels and corresponding face edges

    for idx, polygon in enumerate(polygons):
        center = np.mean(polygon, axis=0)
        dual_node_positions[idx] = tuple((int(center[0]),int(center[1])))
        face_edges = [(node, node + 1) for node in range(len(polygon) - 1)] + [(len(polygon) - 1, 0)]
        original_node_edges = [(res[(int(polygon[edge[0]][0]),int(polygon[edge[0]][1]))], res[(int(polygon[edge[1]][0]),int(polygon[edge[1]][1]))]) for edge in face_edges]
        dual_node_edges[idx] = original_node_edges
        dual_graph.add_node(idx)

    return dual_graph, dual_node_positions, dual_node_edges

# Inside your main code
if __name__ == '__main__':
    polygons = parse_svg('Delta/output-6.svg')

    # Create the dual graph and dictionaries
    dual_graph, d_pos, PL = create_dual_graph(G, polygons)
    print(PL)

    print(d_pos)
    Ge = nx.Graph()
    for u in PL.values():
        Ge.add_edges_from(u)
    Gd = nx.Graph()
    Gd.add_nodes_from(d_pos)
    nx.draw(Ge,pos=Pos,node_size=0)
    nx.draw(Gd,pos=d_pos,node_size=5)
    plt.show()

    # Save the dual graph and dictionaries to files if needed
    pickle.dump(Pos, open('Delta/Pos_output-6_fkt.pickle', 'wb'))
    pickle.dump(d_pos, open('Delta/dual_pos_output-6_fkt.pickle', 'wb'))
    pickle.dump(PL, open('Delta/Plaqs_output-6_fkt.pickle', 'wb'))
'''
# MM = MaxMatch(G,10,Nodes,N*2)
# m = nx.max_weight_matching(G)
# EC = []
# W = []
# for u,v in G.edges():
# 	if((u,v) in m or (v,u) in m):
# 		EC.append('purple')
# 		W.append(2)
# 	else:
# 		EC.append('black')
# 		W.append(0.7)
# nx.draw(G, pos=Pos,node_size=0,edge_color=EC,width=W)   #with_labels=True,font_size=8) #
# #nx.draw(gb, pos=Pos,node_size=10,node_color='orange',with_labels=True,font_size=10)
# plt.show()

'''
c=0
DDL = []
for i in range(10,11):
	M = pickle.load(open('SpectreGamma_MaxMatches_3out/MaxMatch_%i.txt'%i,'rb'))
	l = list(range(90,100))
	#l.remove(i)
	for j in l:
		r=[]
		E = M.copy()
		M1 = pickle.load(open('SpectreGamma_MaxMatches_3out/MaxMatch_%i.txt'%j,'rb'))
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
		nx.draw(G, node_size= 0, pos=Pos, edge_color=EC, width=W) 
		plt.show()
	
#pickle.dump(DDL, open('DDLs_MaxMatchU1_1M_relabel.txt','wb'))
'''

Dim_E = []
EF=[]
for i in range(10):
	E=[]
	MN=[]
	m = pickle.load(open('Delta/MaxMatches_3out_NB/MaxMatch_%i.txt'%i,'rb'))
	for u,v in m:
		MN.append(u)
		MN.append(v)
	DF  = set(Nodes) - set(MN)
	# print(DF)
	Dim_E.append(set(m))
	M_set = set(map(tuple, sorted(m)))
	for u,v in G.edges():
		if ((u,v) not in M_set and (v,u) not in M_set):
			E.append((u,v))
	# for u,v in G.edges():
	# 	if((u,v) not in m and (v,u) not in m):
	# 		E.append((u,v))
	EF.append(set(E))
	# EC = []
	# W1=[]
	# for edge in G.edges():
	# 	if edge in m or tuple(reversed(edge)) in m:   #M_set
	# 		EC.append('purple')
	# 		W1.append(1)
	# 	else:
	# 		EC.append('black')
	# 		W1.append(0.1)
	# pylab.ion()
	# plt.clf()
	# nx.draw(G, pos=Pos, node_size=0, edge_color=EC, width=W1)
	# pause(2)
	# pylab.show()
	#plt.show()

com_EF = set.intersection(*Dim_E)
com_EF1 = set.intersection(*EF)
com = com_EF.copy()
com.update(com_EF1)
Rem = []
for u,v in G.edges():
	if((u,v) not in com and (v,u) not in com):
			Rem.append((u,v))
#pickle.dump(Rem,open('Delta/RedLoops_output-6_fkt.txt','wb'))
#print(len(com_EF),len(com_EF1),len(Rem),len(G.edges()))
EC = []
W1 = []
for u,v in G.edges():
	if((u,v) in com_EF or (v,u) in com_EF):
		EC.append('purple')
		W1.append(1.5)
	elif((u,v) in com_EF1 or (v,u) in com_EF1):
		EC.append('black')
		W1.append(0.7)
	elif((u,v) in Rem or (v,u) in Rem):
		EC.append('orange')
		W1.append(1.5)
	else:
		EC.append('black')
		W1.append(0.7)
	
nx.draw(G, pos=Pos, edge_color=EC, width=W1, node_size=0) #, with_labels=True, font_size=8)
#plt.savefig('X6.svg', format='svg')
#plt.savefig('dimers3.pdf', format='pdf')
plt.show()


Dsgfdgfd

Rem = pickle.load(open('Delta/RedLoops_output-6_fkt.txt','rb'))
g=nx.Graph()
g.add_edges_from(Rem)
nx.draw(g, pos=Pos, node_size=5) #, with_labels=True, font_size=8)
plt.savefig('X6.svg', format='svg')
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

nx.write_gpickle(G,'Spectre_output-5.gpickle')
pickle.dump(Pos, open('Spectre_Pos.txt','wb'))
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
	
