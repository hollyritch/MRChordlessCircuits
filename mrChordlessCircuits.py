from collections import defaultdict
import networkx as nx
from itertools import product
import time
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

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


def chordless_cycle_search_Species(F, B, path, length_bound, G):
    fwBlocked = defaultdict(int)
    bwBlocked = defaultdict(int)
    target = path[0]
    second = path[1]
    for i in range(len(path)):
        a = path[i]
        if i==0 or i==2:
            Fa = F[a]
            for r in Fa:                                            # Block reactions from being visited twice
                fwBlocked[r] += 1
        else:
            Ba=B[a]
            for m in Ba:
                bwBlocked[m] 
    stack = [iter(F[path[2]])]
    while stack:
        nbrs = stack[-1]
        for x in nbrs:
            if x==target:
                continue
            if G.nodes[x]["Type"] == "Species":
                proceed = (bwBlocked[x]==0 and (length_bound is None or len(path) < length_bound))
            else:
                proceed = (fwBlocked[x] == 1 and (length_bound is None or len(path) < length_bound))
            if not proceed:
                continue
            Fx = F[x]
            Bx = B[x]
            if target in Fx:                        # x is a reaction, target a metabolite
                if bwBlocked[target]==0:
                    yield path + [x]                    # yield the path
                    for m in Bx:
                        bwBlocked[m] += 1
                    path.append(x)
                    stack.append(iter(Fx))
                    break
            else:
                if target in Bx:                    # since target is a metabolite, then x is a reaction
                    continue                        # we can probably remove this check, should have been done before
                else:        
                    if G.nodes[x]["Type"]=="Species":
                        if second in Fx and fwBlocked[second]==1:
                            continue
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
            if G.nodes[z]["Type"]=="Species":
                Fz=F[z]
                for y in Fz:
                    fwBlocked[y] -= 1
            else:
                Bz = B[z]
                for y in Bz:
                    bwBlocked[y] -= 1
#############################
#############################


def chordless_cycle_search_Reaction(F, B, path, length_bound, G):
    ## Reaction is the first node
    fwBlocked = defaultdict(int)                                    # will contain only reactions
    bwBlocked = defaultdict(int)                                    # will contain only metabolites‚
    target = path[0]                                                # Target is a reaction
    for i in range(1,len(path)):                                      # Initializing
        a = path[i]
        if i%2==0:                                                  # Reaction case 
            Ba=B[a]
            for m in Ba:
                bwBlocked[m] +=1                                    # Block m from being visited (only unvisited metabolites can be visited)
        else:      
            Fa=F[a]                                                 # Species case
            for r in Fa:                                            # Block reactions from being visited twice
                fwBlocked[r] += 1
    stack = [iter(F[path[2]])]
    while stack:                                                    # Now the real loop
        nbrs = stack[-1]
        for x in nbrs:
            if G.nodes[x]["Type"] == "Species":
                proceed = (bwBlocked[x]==0 and (length_bound is None or len(path) < length_bound))
            else:
                proceed = (fwBlocked[x] == 1 and (length_bound is None or len(path) < length_bound))
            if not proceed:
                continue
            Fx = F[x]
            if target in Fx:                            # Remember: x is a metabolite
                if sum(1 for y in Fx if y in path)==1 and fwBlocked[target]==0:
                    yield path + [x]                    # yield the path and stop.....you can' go any further, otherwise there would be non MR-chordfree cycles
            else:                                   
                if G.nodes[x]["Type"]=="Species":   # Block forward in the case of a metabolite
                    for y in Fx:
                        fwBlocked[y] += 1
                else:                               
                    Bx = B[x]
                    for y in Bx:                    # Block backward in case of a reaction
                        bwBlocked[y] += 1
                path.append(x)                      # Either way, append x to the path add all forward oriented vertices 
                                                    # to the stack
                stack.append(iter(Fx))
                break
        else:                                       # Take off 
            stack.pop() 
            z = path.pop()
            if G.nodes[z]["Type"]=="Species":
                Fz=F[z]
                for y in Fz:
                    fwBlocked[y] -= 1
            else:
                Bz = B[z]
                for y in Bz:
                    bwBlocked[y] -= 1
#############################
#############################


def chordless_cycle_search(F, B, path, length_bound, G):
    if G.nodes[path[0]]["Type"]=="Species":
        yield from chordless_cycle_search_Species(F, B, path, length_bound, G)
    else:
        yield from chordless_cycle_search_Reaction(F, B, path, length_bound, G)
#############################
#############################


def findAllMRChordlessCyclesList(F, bound):
    B = F.reverse(copy=True)
    def stems(C, v):
        for u, w in product(C.pred[v], C.succ[v]):
            yield [u, v, w]

    components = [c for c in nx.strongly_connected_components(F) if len(c) > 3]
    while components:
        c = components.pop()
        v = next(iter(c))   
        Fc = F.subgraph(c)
        Fcc = Bcc = None
        for S in stems(Fc, v):
            if Fcc is None:
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(B.subgraph(c))
            yield from chordless_cycle_search(Fcc, Bcc, S, bound, F)
        components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesSet(F, bound):
    B = F.reverse(copy=True)
    def stems(C, v):
        for u, w in product(C.pred[v], C.succ[v]):
            yield [u, v, w]

    components = [c for c in nx.strongly_connected_components(F) if len(c) > 3]
    while components:
        c = components.pop()
        v = next(iter(c))   
        Fc = F.subgraph(c)
        Fcc = Bcc = None
        for S in stems(Fc, v):
            if Fcc is None:
                Fcc = _NeighborhoodCache(Fc)
                Bcc = _NeighborhoodCache(B.subgraph(c))
            yield from chordless_cycle_search(Fcc, Bcc, S, bound, F)
        components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
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


def callFindMRChordlessCircuits(H:nx.DiGraph, type:str, bound:int):
    for comp in nx.strongly_connected_components(H):
        F = H.subgraph(comp)
        X, Y = nx.algorithms.bipartite.sets(F, top_nodes=None)
        if len(X)<len(Y):
            M=X
            R=Y
        else:
            M=Y
            R=X
        for m in M:
            F.nodes[m]["Type"]="Species"
        for r in R:
            F.nodes[r]["Type"]="Reaction"
        for n in F.nodes():
            if "Type" not in F.nodes[n]:
                print(n, n in X, n in R, F.nodes[n])
                input()
        if type == "List":
            for c in findAllMRChordlessCyclesList(F, bound):
                continue
        if type == "Set":
            for c in findAllMRChordlessCyclesSet(F, bound):
                continue
#############################
#############################


def checkCycle(c, H):
    chordFree = True
    for v in c:
        if H.nodes[v]["Type"]=="Reaction":            
            if sum(1 for u in c if (u,v) in H.edges())>1:
                chordFree=False
                break
        if H.nodes[v]["Type"]=="Species":            
            if sum(1 for u in c if (v,u) in H.edges())>1:
                chordFree=False
                break
    return chordFree
#############################
#############################


def plotResults(m, p, keys, timeDict):
    x = np.arange(len(keys))  # the label locations
    width = 0.1  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')
    fig.set_figheight(10)
    fig.set_figwidth(20)
    for attribute, measurement in timeDict.items():
        offset = width * multiplier
        rects = ax.bar(x+0.1 + offset, measurement, width, label=attribute)
        multiplier += 1
    n = 28
    ax.set_ylabel(ylabel="Time", fontsize=n)
    ax.set_yscale("log")
    ax.set_xlabel(xlabel="Number of nodes", fontsize=n)
    ax.set_xticks(x + width, keys, fontsize=n)
    ax.legend(loc='upper left', ncols=3, fontsize=n)
    plt.yticks(fontsize=n)
    plt.savefig("./Benchmarking/MRChordlessVsJohnson_"+str(p)+".png")
#############################
#############################

x = 10
m = 10
bound=20
for k in tqdm(range(1, 101, 1), leave=False, total=20):
    p = k/100
    print(p)
    timeDict = {}
    keys = []
    list = range(m)
    for i in tqdm(range(len(list)), leave = False, total=len(list)):
        if i==0:
            continue
        n=i*x
        H=nx.algorithms.bipartite.random_graph(n=n, m=n, p=p, directed=True)
        
        timeStamp = time.time()
        print("Starting to compute MR-Chordless cycles with set")
        callFindMRChordlessCircuits(H, "Set", bound)
        setTime = time.time()-timeStamp
        
        timeStamp = time.time()
        print("Starting to compute MR-Chordless cycles with list")
        callFindMRChordlessCircuits(H, "List", bound)
        listTime = time.time()-timeStamp
        
        timeStamp = time.time()
        print("Starting to compute all cycles")
        for c in nx.simple_cycles(H, length_bound=bound):
            checkCycle(c, H)
        johnsonTime = time.time()-timeStamp
        if n==1*x:
            timeDict["MR-chordlessSet"]=[setTime]
            timeDict["MR-chordlessList"]=[listTime]
            timeDict["Johnson"]=[johnsonTime]
        else:
            timeDict["MR-chordlessSet"].append(setTime)
            timeDict["MR-chordlessList"].append(listTime)
            timeDict["Johnson"].append(johnsonTime)
        keys.append(n)
    plotResults(2*n, p, keys, timeDict)