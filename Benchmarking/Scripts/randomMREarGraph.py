#!/usr/bin/env python3
import networkx as nx
#from pyvis.network import Network
from random import choice
from random import randint


def add_ear(index, length, G, n1, n2, X, R):
    next_is_R = G.nodes[n1]['M']   # if starting node is a compound, start with reaction
    last_node = n1    
    for i in range(length):
        new_node = "E"+str(index)+"_"+str(len(G.nodes()))        
        G.add_node(new_node, R=next_is_R, M = not next_is_R)
        if G.nodes[new_node]["M"]==True:
            X.add(new_node)
        else:
            R.add(new_node)
        G.add_edge(last_node, new_node)
        next_is_R = not next_is_R
        last_node = new_node
        
    if (G.nodes[n1]['M'] == G.nodes[n2]['M'] and length %2==0) or (G.nodes[n1]['M'] != G.nodes[n2]['M'] and length%2!=0): # add one more to keep bipartite

        new_node = "E"+str(index)+"_"+str(len(G.nodes()))
        G.add_node(new_node, R=next_is_R, M=not next_is_R)
        if G.nodes[new_node]["M"]==True:
            X.add(new_node)
        else:
            R.add(new_node)
        G.add_edge(last_node, new_node)
        next_is_R = not next_is_R        
        last_node = new_node   
        
    G.add_edge(last_node, n2)
    if G.nodes[n1]["R"]==False or G.nodes[n2]["M"]==False:
        print("Somethings wrong")
        input()

#############################
#############################    


def split_edge(G, edge, X, R, start):
    n1 = edge[0]  
    n2 = edge[1]
  
    new1 = "S"+str(len(G.nodes()))
    new2 = "S"+str(len(G.nodes())+1)

    G.add_node(new1, R=G.nodes[n1]['M'], M=G.nodes[n1]['R'])
    G.add_node(new2, R=G.nodes[n1]['R'], M=G.nodes[n1]['M'])
    
    G.add_edge(n1, new1)
    G.add_edge(new1, new2)
    G.add_edge(new2, n2)
    
    if start==True:
        if G.nodes[n1]["R"]==True:
            X.add(new1)
            R.add(new2)
            return n1
        else:
            X.add(new2)
            R.add(new1)
            return n2
    else:
        if G.nodes[n1]["R"]==True:
            X.add(new1)
            R.add(new2)
            return n2
        else:
            X.add(new2)
            R.add(new1)
            return n1

    #return choice([n1,n2])
#############################
#############################    


# Define initial graph
def mainRandomMREarGraph(ears_to_add, max_ear_length):
    G = nx.DiGraph()

    G.add_node("R_1", R=True, M=False)
    G.add_node("R_2", R=True, M=False)
    G.add_node("M_1", R=False, M=True)
    G.add_node("M_2", R=False, M=True)
    
    X = {n for n,d in G.nodes(data=True) if d["M"]==True}
    R = set(G)-X

    G.add_edge("M_1", "R_1")
    G.add_edge("R_1", "M_2")
    G.add_edge("M_2", "R_2")
    G.add_edge("R_2", "M_1")
    # Add ears
    for i in range(ears_to_add):
        
        start = choice( list(R)+ list(G.edges))
        end = choice( list(X) + list(G.edges))
        
        startnode = start 
        if start not in G.nodes():
            while G.nodes[start[0]]["R"]==False:
                start=choice(list(G.edges))
            startnode = split_edge(G, start, X, R, True)                
            
            endnode = end 
            if end not in G.nodes():
                while G.nodes[end[0]]["R"]==False:
                    end=choice(list(G.edges))
                endnode = split_edge(G, end, X, R, False)            
        else:   
            endnode = end
            if end not in G.nodes():
                while G.nodes[end[0]]["R"]==False:
                    end=choice(list(G.edges))
                endnode = split_edge(G, end, X, R, False)
        add_ear(i, randint(2, max_ear_length), G, startnode, endnode, X, R) 
    return G, X, R
#############################
#############################    

