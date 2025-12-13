from collections import defaultdict
import networkx as nx
from itertools import product
import time
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
import pickle
import sys, os
from copy import deepcopy

class _NeighborhoodCache(dict):
    """Very lightweight graph wrapper which caches neighborhoods as list.

    This dict subclass uses the __missing__ functionality to query graphs for
    their neighborhoods, and store the result as a list.  This is used to avoid
    the performance penalty incurred by subgraph views.
    """

    def __init__(self, G):
        self.G = G

    def __missing__(self, v):
        Gv = self[v] = set(self.G[v])
        return Gv
#############################
#############################


def chordless_cycle_search_Species(F:nx.DiGraph, B:nx.DiGraph, path:list, length_bound:int, metabolites:set):
    fwBlocked = defaultdict(int)
    bwBlocked = defaultdict(int)
    target = path[0]
    for i in range(len(path)):
        a = path[i]
        if i==0 or i==2:
            Fa = F[a]
            for r in Fa:                                            # Block reactions from being visited twice
                fwBlocked[r] += 1
        else:
            Ba=B[a]
            for m in Ba:
                bwBlocked[m] +=1

    stack = [iter(F[path[2]])]
    if fwBlocked[path[1]]>1:
        return
    while stack:
        nbrs = stack[-1]
        for x in nbrs:
            if x==target:
                continue
            if x in metabolites:
                proceed = (bwBlocked[x]==0 and (length_bound is None or len(path) < length_bound))
            else:
                proceed = (fwBlocked[x] == 1 and (length_bound is None or len(path) < length_bound))
            if not proceed:
                continue
            Fx = F[x]
            if target in Fx:                        # x is a reaction, target a metabolite    
                if bwBlocked[target]==1:                    
                    yield path + [x]
                    Bx = B[x]
                    for m in Bx:
                        bwBlocked[m] += 1
                    path.append(x)
                    stack.append(iter(Fx))
                    break
            else:
                Bx = B[x]
                if target in Bx:                    # since target is a metabolite, then x is a reaction
                    continue                        # we can probably remove this check, should have been done before
                else:        
                    if x in metabolites:
                        for r in Fx:
                            fwBlocked[r] += 1
                    else:
                        for m in Bx:
                            bwBlocked[m] += 1
                    path.append(x)
                    stack.append(iter(Fx))
                    break
        else:                                       # Take off 
            stack.pop() 
            z = path.pop()
            if z in metabolites:
                Fz=F[z]
                for y in Fz:
                    fwBlocked[y] -= 1
            else:
                Bz = B[z]
                for y in Bz:
                    bwBlocked[y] -= 1
#############################
#############################


def findAllMRChordlessCyclesList(F, reactions:set, metabolites:set, bound:int):
    B = F.reverse(copy=True)
    for u in metabolites:
        Fu = F.successors(u)
        digons = [[u, v] for v in Fu if F.has_edge(v, u)]
        yield from digons
        
    def stems(C, v):
        for u, w in product(C.pred[v], C.succ[v]):
            yield [u, v, w]

    components = [c for c in nx.strongly_connected_components(F) if len(c) > 3]
    while components:
        c = components.pop()
        rx = reactions & c
        if rx:
            v = next(iter(rx))
            Fc = F.subgraph(c)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesSet(F, reactions:set, metabolites:set, bound:int):
    B = F.reverse(copy=True)
    for u in metabolites:
        Fu = F.successors(u)
        digons = [[u, v] for v in Fu if F.has_edge(v, u)]
        yield from digons
        
    def stems(C, v):
        for u, w in product(C.pred[v], C.succ[v]):
            yield [u, v, w]

    separate = nx.strongly_connected_components
    components = [c for c in separate(F) if len(c) > 3]
    while components:
        c = components.pop()
        rx = reactions & c
        if rx:
            v = next(iter(rx))
            Fc = F.subgraph(c)
            Fcc = Bcc = None
            for S in stems(Fc, v):
                if Fcc is None:
                    Fcc = _NeighborhoodCache(Fc)
                    Bcc = _NeighborhoodCache(B.subgraph(c))
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in separate(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findNodeSpecificMRChordlessCycles(F, x, bound):
    B = F.reverse()
    def stems(C, v):
        for u, w in product(C.pred[v], C.succ[v]):
            yield [u, v, w]

    components = [c for c in nx.strongly_connected_components(F) if len(c) > 3]
    while components:
        c = components.pop()
        if x in c:
            Fc = F.subgraph(c)
            Bc = B.subgraph(c)
            Fcc = Bcc = None
            for S in stems(Fc, x):
                if Fcc is None:
                    Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                    Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search(Fcc, Bcc, S, bound, G)
            break
#############################
#############################


def checkCycle(c, K):
    chordFree = True
    for v in c:
        if K.nodes[v]["Type"]=="Reaction":            
            if sum(1 for u in c if H.has_edge(u,v))>1:
                chordFree=False
                break
        if K.nodes[v]["Type"]=="Species":            
            if sum(1 for u in c if H.has_edge(v,u))>1:
                chordFree=False
                break
    return chordFree
#############################
#############################


def plotResults(minNoNodes:int, maxNoNodes:int, p:float, keys:set, path:str, timeDict:dict):
    x = np.arange(len(keys))  # the label locations
    width = 0.1  # the width of the bars
    multiplier = 0
    newTimeDict = {}
    i=0
    for counter in keys:
        setList= timeDict["MR-chordlessSet"][counter]
        listList= timeDict["MR-chordlessList"][counter]
        johnsonList= timeDict["Johnson"][counter]
        newTimeDict.setdefault("MR-chordlessSet", []).append(sum(setList)/len(setList))
        newTimeDict.setdefault("MR-chordlessList", []).append(sum(listList)/len(listList))
        newTimeDict.setdefault("Johnson", []).append(sum(johnsonList)/len(johnsonList))
    
    fig, ax = plt.subplots(layout='constrained')
    fig.set_figheight(10)
    fig.set_figwidth(20)
    for attribute, measurement in newTimeDict.items():
        offset = width * multiplier
        rects = ax.bar(x+0.1 + offset, measurement, width, label=attribute)
        multiplier += 1
    n = 28
    ax.set_ylabel(ylabel="Time", fontsize=n)
    ax.set_title(label="Random graphs of max. size " +str(maxNoNodes), fontsize=n+4)
    ax.set_yscale("log")
    ax.set_xlabel(xlabel="Number of MR-chordless cycles", fontsize=n)
    ax.set_xticks(x + width, keys, fontsize=n, rotation=45)
    ax.legend(loc='upper left', ncols=1, fontsize=n)
    plt.yticks(fontsize=n)
    plt.savefig(path + "/MRChordlessVsJohnson_"+str(p)+"min_" + str(minNoNodes) + "max_" +str(maxNoNodes)+".png")
    plt.close()
#############################
#############################


def buildExampleGraph():
    G=nx.DiGraph()
    X=set()
    R=set()

    # Metabolites
    G.add_node(7)
    G.add_node(27)
    G.add_node(14)
    G.add_node(19)
    X.add(7)
    X.add(27)
    X.add(14)
    X.add(19)

    # Reactions
    G.add_node(36)
    G.add_node(37)
    G.add_node(58)
    G.add_node(55)
    R.add(36)
    R.add(37)
    R.add(58)
    R.add(55)

    G.add_edge(7,58)
    G.add_edge(7,36)
    G.add_edge(36,27)
    G.add_edge(27,37)
    G.add_edge(37,14)
    G.add_edge(14,58)
    G.add_edge(58,19)
    G.add_edge(19,55)
    G.add_edge(55,7)
    return G, X, R

l = int(sys.argv[1])
maxSize = 6

a = float(sys.argv[2])
b = float(sys.argv[3])

direc= "./Benchmarking/"+str(a)+"to"+str(b)
if not os.path.exists(direc):
    os.makedirs(direc)

bound=None
timeDict = {}
cycleDict = {}
myList = range(maxSize)

G,X,R = buildExampleGraph()


for k in tqdm(range(1, 11, 1), leave=False, total=10):
    p = k/200
    keys = set()
    for i in tqdm(range(1,len(myList)), leave = False, total=len(myList)):
        noMetabolites=int(a*i*l)
        noReactions=int(b*i*l)
        if i ==1:
            minNoNodes = noMetabolites+noReactions
        elif i==len(myList)-1:
            maxNoNodes = noMetabolites+noReactions
        for j in tqdm(range(l), leave = False, total=l):
            G=nx.algorithms.bipartite.random_graph(n=noMetabolites, m=noReactions, p=p, directed=True)
            X = {n for n,d in G.nodes(data=True) if d['bipartite'] == 0}
            R = set(G) - X
            metabolites = set()
            reactions = set()
            # assign labels
            for x in X:
                G.nodes[x]['Type'] = "Species"
            for r in R:
                G.nodes[r]['Type'] = "Reaction"
            
            H=deepcopy(G)
            K=deepcopy(G)
            L=deepcopy(G)

            timeStamp = time.time()
            setCounter = sum(1 for c in findAllMRChordlessCyclesSet(H, R, X, bound))
            setTime = time.time()-timeStamp
            
            timeStamp = time.time()
            listCounter = sum(1 for c in findAllMRChordlessCyclesList(K, R, X, bound))
            listTime = time.time()-timeStamp
            
            if setCounter!=listCounter:
                print("Set and list counter are not the same", setCounter, listCounter)
                for c in findAllMRChordlessCyclesSet(deepcopy(G), R, X, bound):
                    print(c)
                print()
                for c in findAllMRChordlessCyclesList(deepcopy(G), R, X, bound):
                    print(c)
                input()

            timeStamp = time.time()
            johnsonCounter = 0
            totalJohnsonCounter = 0
            for c in nx.simple_cycles(L, length_bound=bound):
                totalJohnsonCounter+=1
                if checkCycle(c, L):
                    johnsonCounter+=1
            if setCounter!=johnsonCounter:                
                print("Set and Johnson counter are not the same", setCounter, johnsonCounter)
                print("List and Johnson-counter are", listCounter, johnsonCounter)
                for c in findAllMRChordlessCyclesList(deepcopy(G), R, X, bound):
                    print(c)
                print()
                for c in nx.simple_cycles(K, length_bound=bound):
                    if checkCycle(c, K):
                        print(c)

                x = int(input())
                for c in findAllMRChordlessCyclesList2(deepcopy(G), R, X, bound, x):
                    print("New cycle is", c)
            johnsonTime = time.time()-timeStamp
            if "MR-chordlessSet" not in timeDict:
                timeDict["MR-chordlessSet"]={setCounter: [setTime]}
            if "MR-chordlessList" not in timeDict:
                timeDict["MR-chordlessList"]={setCounter: [listTime]}
            if "Johnson" not in timeDict:
                timeDict["Johnson"]={setCounter: [johnsonTime]}
            else:
                timeDict["MR-chordlessSet"].setdefault(setCounter, []).append(setTime)
                timeDict["MR-chordlessList"].setdefault(setCounter, []).append(listTime)
                timeDict["Johnson"].setdefault(setCounter, []).append(johnsonTime)
            keys.add(setCounter)
            cycleDict[setCounter]=totalJohnsonCounter
    path = direc+ "/Benchmarking" +str(p) + "min_" + str(minNoNodes) + "max_" +str(maxNoNodes) +".pkl"
    plotResults(minNoNodes, maxNoNodes, p, sorted(list(keys)), direc, timeDict)
    with open(path, "wb") as file:
        pickle.dump((timeDict, cycleDict), file)
