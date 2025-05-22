import matplotlib.pylab as plt
import networkx as nx
import numpy as np
import random
import pickle

class ASI():
    def angle(self,e1,e2):
        a = np.array(pos[e1[0]]) - np.array(pos[e1[1]])
        b = np.array(pos[e2[0]]) - np.array(pos[e2[1]])
        cosi = round(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)),1)
        if cosi == 1 or cosi == -1:
            return True
        else:
            return False
    def directions(self):
        DirE = []
        node = random.choice(list(G.nodes()))
        d = 1
        Vn = []
        flag = False
        while len(DirE) != len(G.edges()):
            #print(len(DirE))
            #print(node)
            if G.degree(node) == 2:
                for i in G.neighbors(node):
                    if (i,node) not in DirE and (node,i) not in DirE:
                        flag = True
                        if d == 1:
                            #print('*')
                            DirE.append((node,i))
                        else:
                            #print('**')
                            DirE.append((i,node))
                        node = i
                        break
            else:
                for i in G.neighbors(node):
                    if (i,node) not in DirE and (node,i) not in DirE:
                        Vn.append(node)
                        flag = True
                        if DirE == []:
                            DirE.append((node,i))
                        elif self.angle(pe,(node,i)):
                            if d == 1:
                                #print('#')
                                d = -1
                                DirE.append((i,node))
                            else:
                                #print('##')
                                d = 1
                                DirE.append((node,i))
                        else:
                            if d == 1:
                                #print('%')
                                DirE.append((node,i))
                            else:
                                #print('%%')
                                DirE.append((i,node))
                        node = i
                        break
            if flag == False:
                node = random.choice(Vn)
                #print('nn',node)
                for j in DirE:
                    if node in j:
                        pe = j
                        break
                #print(pe)
                for i in G.neighbors(node):
                    if (i, node) not in DirE and (node, i) not in DirE:
                        Vn.append(node)
                        flag = True
                        if self.angle(pe, (node, i)):
                            if pe[1] == node:
                                #print('&')
                                d = -1
                                DirE.append((i, node))
                            else:
                                #print('&&')
                                d = 1
                                DirE.append((node, i))
                        else:
                            if pe[1] == node:
                                #print('£')
                                d = 1
                                DirE.append((node, i))
                            else:
                                #print('££')
                                d = -1
                                DirE.append((i, node))
                        node = i
                        break
            else:
                pe = DirE[-1]
                flag = False
        return DirE

G = pickle.load(open('Hout2.txt', 'rb'))
pos = pickle.load(open('Pos_Hout2.txt', 'rb'))

obj = ASI()
DirE = obj.directions()

g = nx.DiGraph()
g.add_edges_from(DirE)
print(DirE)
nx.draw(g,pos,node_size=5)  #,with_labels=True)
plt.show()

