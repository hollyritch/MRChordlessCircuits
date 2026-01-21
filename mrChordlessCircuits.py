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
import randomMREarGraph
from itertools import cycle
from itertools import chain


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


def chordless_cycle_search_Species(F:nx.DiGraph, B:nx.DiGraph, path:list, length_bound:int, X:set):
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
            if x in X:
                proceed = (bwBlocked[x]==0 and (length_bound is None or len(path) < length_bound))
            else:
                proceed = (fwBlocked[x] == 1 and (length_bound is None or len(path) < length_bound))
            if not proceed:
                continue
            Fx = F[x]
            if target in Fx:                        # x is a reaction, target a metabolite    
                #if bwBlocked[target]==1:                    
                yield path + [x]
                Bx = B[x]
                for m in Bx:
                    bwBlocked[m] += 1
                path.append(x)
                stack.append(iter(Fx))
                break
            else:
                Bx = B[x]
                # if target in Bx:                    # since target is a metabolite, then x is a reaction
                #     continue                        # we can probably remove this check, should have been done before
                # else:        
                if x in X:
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
            if z in X:
                Fz=F[z]
                for y in Fz:
                    fwBlocked[y] -= 1
            else:
                Bz = B[z]
                for y in Bz:
                    bwBlocked[y] -= 1
#############################
#############################


def findAllMRChordlessCyclesListSimpleDegree(F, reactions:set, metabolites:set, bound:int):
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
            Fc = F.subgraph(c)
            v = max(rx, key = Fc.degree)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesListInDegree(F, reactions:set, metabolites:set, bound:int):
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
            Fc = F.subgraph(c)
            v = max(rx, key = Fc.in_degree)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesListOutDegree(F, reactions:set, metabolites:set, bound:int):
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
            Fc = F.subgraph(c)
            v = max(rx, key = Fc.out_degree)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesListReacInDegree(F, reactions:set, metabolites:set, bound:int):
    B = F.reverse(copy=True)
    R = createReactionNetwork(F, reactions)
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
            Rc = R.subgraph(rx)
            v = max(rx, key = Rc.in_degree)
            Fc = F.subgraph(c)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesListReacOutDegree(F, reactions:set, metabolites:set, bound:int):
    B = F.reverse(copy=True)
    R = createReactionNetwork(F, reactions)
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
            Rc = R.subgraph(rx)
            v = max(rx, key = Rc.out_degree)
            Fc = F.subgraph(c)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findAllMRChordlessCyclesListReacSymmDegree(F, reactions:set, metabolites:set, bound:int):
    B = F.reverse(copy=True)
    R = createReactionNetwork(F, reactions)
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
            Rc = R.subgraph(rx)
            v = max(rx, key = Rc.degree)
            Fc = F.subgraph(c)
            Bc = B.subgraph(c)
            for S in stems(Fc, v):
                if S[0]==S[2]:
                    continue
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
                yield from chordless_cycle_search_Species(Fcc, Bcc, S, bound, metabolites)
            components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
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
                if S[0]==S[2]:
                    continue
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


def createReactionNetwork(G:nx.DiGraph, reactions:set):
    reactionNetwork = nx.DiGraph()
    for r1 in reactions:
        for r2 in reactions:
            if r1==r2:
                continue
            # First Check if a product of r1 is a reactant of r2
            productsR1 = set(G.successors(r1))
            eductsR2 = set(G.predecessors(r2))
            if len(productsR1.intersection(eductsR2))>0:
                reactionNetwork.add_edge(r1, r2)
    return reactionNetwork
#############################
#############################


def checkCycle(c, G, X, R):
    chordFree = True
    for v in c:
        if v in R:            
            if sum(1 for u in c if G.has_edge(u,v))>1:
                chordFree=False
                break
        if v in X:            
            if sum(1 for u in c if G.has_edge(v,u))>1:
                chordFree=False
                break
    return chordFree
#############################
#############################


def plotResultsBar(minNoNodes:int, maxNoNodes:int, p:float, keys:set, path:str, timeDict:dict):
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
    plt.savefig(path + "/MRChordlessVsJohnson_Bar_"+str(p)+"min_" + str(minNoNodes) + "max_" +str(maxNoNodes)+".png")
    plt.close()
#############################
#############################


def plotResultsLine(minNoNodes:int, maxNoNodes:int, p:float, path:str, timeDict:dict):
    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(30)
    barLabel=["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink", "tab:brown", "tab:gray"]
    i = 0
    n=32
    for algo, sizeDict in timeDict.items():
        keys = sorted(list(sizeDict.keys()))
        values = []
        col = barLabel[i]
        for l in keys:
            values.append(sum(sizeDict[l])/len(sizeDict[l]))
        i+=1
        ax.plot(keys, values, linewidth=1, label=algo, color=col)
    ax.legend(loc='upper left', ncols=1, fontsize=n)
    ax.set_title("Enumeration of MR chordfree elementary circuits vs. Johnson algorithm", fontsize=n+4)
    ax.set_xlabel("Number of MR chordfree circuits", fontsize=n)
    ax.set_ylabel("Time in seconds", fontsize=n)
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    plt.savefig(path + "/MRChordlessVsJohnson_Line_"+str(p)+"min_" + str(minNoNodes) + "max_" +str(maxNoNodes)+".png")
    plt.close()
    return
#############################
#############################


def plotResultsLineEars(k:int, length:int, path:str, timeDict:dict):
    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(30)
    barLabel=["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink", "tab:brown", "tab:gray"]
    i = 0
    n=32
    for algo, sizeDict in timeDict.items():
        keys = sorted(list(sizeDict.keys()))
        values = []
        col = barLabel[i]
        for l in keys:
            values.append(sum(sizeDict[l])/len(sizeDict[l]))
        i+=1
        ax.plot(keys, values, linewidth=1, label=algo, color=col)
    ax.legend(loc='upper left', ncols=1, fontsize=n)
    ax.set_xlabel("Number of MR chordfree circuits", fontsize=n)
    ax.set_ylabel("Time in seconds", fontsize=n)
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    plt.savefig(path + "/MRChordlessVsJohnsonVsFluffles_Line_" + str(k) + "_MaxEarLength_" + str(length)+".png")
    plt.close()
    return
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
#############################
#############################


def checkIfCountersMatch(setCounter, listCounter, listDegreeCounter, listROutDegreeCounter, listRSymmDegreeCounter,johnsonCounter):
    if setCounter!=listCounter:
        print("Set and list counter are not the same", setCounter, listCounter)
        for c in findAllMRChordlessCyclesSet(deepcopy(G), R, X, bound):
            print(c)
        print()
        for c in findAllMRChordlessCyclesList(deepcopy(G), R, X, bound):
            print(c)
        input()
    if setCounter != listDegreeCounter:
        print("Set and listDegree counter are not the same", setCounter, listCounter)
        for c in findAllMRChordlessCyclesSet(deepcopy(G), R, X, bound):
            print(c)
        print()
        for c in findAllMRChordlessCyclesListSimpleDegree(deepcopy(G), R, X, bound):
            print(c)
        input()
    if setCounter != listROutDegreeCounter:
        print("Set and listDegree counter are not the same", setCounter, listCounter)
        for c in findAllMRChordlessCyclesSet(deepcopy(G), R, X, bound):
            print(c)
        print()
        for c in findAllMRChordlessCyclesListReacOutDegree(deepcopy(G), R, X, bound):
            print(c)
        input()
    if setCounter != listRSymmDegreeCounter:
        print("Set and listDegree counter are not the same", setCounter, listCounter)
        for c in findAllMRChordlessCyclesSet(deepcopy(G), R, X, bound):
            print(c)
        print()
        for c in findAllMRChordlessCyclesListReacSymmDegree(deepcopy(G), R, X, bound):
            print(c)
        input()
    if setCounter!=johnsonCounter:                
        print("Set and Johnson counter are not the same", setCounter, johnsonCounter)
        print("List and Johnson-counter are", listCounter, johnsonCounter)
        for c in findAllMRChordlessCyclesList(deepcopy(G), R, X, bound):
            print(c)
        print()
        for c in nx.simple_cycles(K, length_bound=bound):
            if checkCycle(c, K):
                print(c)
        input()
    return
#############################
#############################


def assignToTimeDict(timeDict:dict, setCounter:int, setTime:float, listTime:float, listDegreeTime:float, listInDegreeTime:float, listOutDegreeTime:float, listRInDegreeTime:float, listROutDegreeTime:float, listRSymmDegreeTime:float, johnsonTime:float):

    if "MR-chordlessSet" not in timeDict.keys():
        timeDict["MR-chordlessSet"]={}    
    if setCounter in timeDict["MR-chordlessSet"].keys():
        timeDict["MR-chordlessSet"][setCounter].append(setTime)
    else:
        timeDict["MR-chordlessSet"][setCounter]=[setTime]

    if "MR-chordlessList" not in timeDict:
        timeDict["MR-chordlessList"]={}    
    if setCounter in timeDict["MR-chordlessList"].keys():
        timeDict["MR-chordlessList"][setCounter].append(listTime)
    else:
        timeDict["MR-chordlessList"][setCounter]=[listTime]

    if "MR-ChordlessListSimpleDegree" not in timeDict:
        timeDict["MR-ChordlessListSimpleDegree"]={}    
    if setCounter in timeDict["MR-ChordlessListSimpleDegree"].keys():
        timeDict["MR-ChordlessListSimpleDegree"][setCounter].append(listDegreeTime)
    else:
        timeDict["MR-ChordlessListSimpleDegree"][setCounter]=[listDegreeTime]

    if "MR-ChordlessListInDegree" not in timeDict:
        timeDict["MR-ChordlessListInDegree"]={}
    if setCounter in timeDict["MR-ChordlessListInDegree"].keys():
        timeDict["MR-ChordlessListInDegree"][setCounter].append(listInDegreeTime)
    else:
        timeDict["MR-ChordlessListInDegree"][setCounter]=[listInDegreeTime]

    if "MR-ChordlessListOutDegree" not in timeDict:
        timeDict["MR-ChordlessListOutDegree"]={}
    if setCounter in timeDict["MR-ChordlessListOutDegree"].keys():
        timeDict["MR-ChordlessListOutDegree"][setCounter].append(listOutDegreeTime)
    else:
        timeDict["MR-ChordlessListOutDegree"][setCounter]=[listOutDegreeTime]

    if "MR-ChordlessListRInDegreeTime" not in timeDict:
        timeDict["MR-ChordlessListRInDegreeTime"]={}    
    if setCounter in timeDict["MR-ChordlessListRInDegreeTime"].keys():
        timeDict["MR-ChordlessListRInDegreeTime"][setCounter].append(listRInDegreeTime)
    else:
        timeDict["MR-ChordlessListRInDegreeTime"][setCounter]=[listRInDegreeTime]

    if "MR-ChordlessListROutDegreeTime" not in timeDict:
        timeDict["MR-ChordlessListROutDegreeTime"]={}    
    if setCounter in timeDict["MR-ChordlessListROutDegreeTime"].keys():
        timeDict["MR-ChordlessListROutDegreeTime"][setCounter].append(listROutDegreeTime)
    else:
        timeDict["MR-ChordlessListROutDegreeTime"][setCounter]=[listROutDegreeTime]

    if "MR-ChordlessListRSymmDegreeTime" not in timeDict:
        timeDict["MR-ChordlessListRSymmDegreeTime"]={}    
    if setCounter in timeDict["MR-ChordlessListRSymmDegreeTime"].keys():
        timeDict["MR-ChordlessListRSymmDegreeTime"][setCounter].append(listRSymmDegreeTime)
    else:
        timeDict["MR-ChordlessListRSymmDegreeTime"][setCounter]=[listRSymmDegreeTime]

    if "Johnson" not in timeDict:
        timeDict["Johnson"]={}
    if setCounter in timeDict["Johnson"].keys():
        timeDict["Johnson"][setCounter].append(johnsonTime)
    else:
        timeDict["Johnson"][setCounter]=[johnsonTime]
#############################
#############################


def assignToTimeDictEears(timeDict:dict, listCounter:int, listTime:float, johnsonTime:float, fluffleTime:float):
    if "MR-chordlessList" not in timeDict:
        timeDict["MR-chordlessList"]={}
    if listCounter in timeDict["MR-chordlessList"].keys():
        timeDict["MR-chordlessList"][listCounter].append(listTime)
    else:
        timeDict["MR-chordlessList"][listCounter]=[listTime]

    if "FluffleTime" not in timeDict:
        timeDict["FluffleTime"]={}
    if listCounter in timeDict["FluffleTime"].keys():
        timeDict["FluffleTime"][listCounter].append(fluffleTime)
    else:
        timeDict["FluffleTime"][listCounter]=[fluffleTime]
    
    if "Johnson" not in timeDict:
        timeDict["Johnson"]={}
    if listCounter in timeDict["Johnson"].keys():
        timeDict["Johnson"][listCounter].append(johnsonTime)
    else:
        timeDict["Johnson"][listCounter]=[johnsonTime]
#############################
#############################


def assembleFluffles(queue:list, G:nx.DiGraph, X:set, R:set):
    if len(queue)==0:
        return 0
    # 0. Define new variables
    additionalCycleCounter = 0
    additionalCycles = []
    # Generate important dictionaries for further analysis
    edgeCycleDict, cycleIDDict, cycleIDEdgeDict, cycleIDNodeDict, visitedEdges = generateEdgeCycleDict(queue)
    newQueue = deepcopy(list(cycleIDDict.keys()))
    initialNoCycles = len(cycleIDDict.keys())
    spinner = cycle("|/-\\")
    while True:
        cKey = newQueue.pop(0)
        c = cycleIDDict[cKey]
        intersectingCycles, edgesC = getIntersectingCycles(cKey, cycleIDEdgeDict, edgeCycleDict)
        for cIntID in intersectingCycles:
            cInt = cycleIDDict[cIntID]
            if cInt == c:                                                                           # if the cycle is identical to the one we are currently looking at
                continue
            else:
                edgesC2 = cycleIDEdgeDict[cIntID]                                                   # Read edges for the cycle we are looking at 
                allEdges = frozenset(edgesC.union(edgesC2))                                         # Generate set of all edges of c1 and c2                
                if allEdges in visitedEdges:                                                        # If we have seen this set of edges already, then continue
                    continue
                if cKey >= initialNoCycles:
                    cNew = [*c,cInt]
                else:
                    cNew = [c,cInt]
                additionalCycles.append(cNew)
                visitedEdges.add(allEdges)
                plausible = checkPlausibilityOfMatching(G, cKey, cIntID, cycleIDEdgeDict, cycleIDNodeDict, X, R)
                if plausible == True:                    
                    fluffle = True
                    cNewFlat = list(chain.from_iterable(cNew))
                    F = G.subgraph(cNewFlat)
                    for v in cNewFlat:
                        if v in X:
                            if F.out_degree(v)>1:
                                fluffle=False
                                break
                        else:
                            if F.in_degree(v)>1:
                                fluffle=False
                                break
                    if fluffle==True:
                        additionalCycleCounter +=1                                                          # Increase number of analysed cycles
                        newKey = len(cycleIDDict.keys())
                        newQueue.append(newKey)                    
                        cycleIDDict[newKey]= cNew
                        cycleIDEdgeDict[newKey] = cycleIDEdgeDict[cKey].union(cycleIDEdgeDict[cIntID])
                        cycleIDNodeDict[newKey] = cycleIDNodeDict[cKey].union(cycleIDNodeDict[cIntID])
                    else:
                        continue
                else:
                    continue
                    #print("Unfortnuately not plausible")
        sys.stdout.write(f"\r{next(spinner)} Queue length: {len(newQueue):<5}")
        sys.stdout.flush()
        if len(newQueue) == 0:
            break
    return additionalCycleCounter
#############################
#############################  


def checkPlausibilityOfMatching(G:nx.DiGraph, c1Key:int, c2Key:int, cycleIDEdgeDict:dict, cycleIDNodeDict:dict, X:set, R:set):
    edgesC1 = cycleIDEdgeDict[c1Key]
    edgesC2 = cycleIDEdgeDict[c2Key]
    nodesC1 = cycleIDNodeDict[c1Key]
    nodesC2 = cycleIDNodeDict[c2Key]
    intersectingEdges = edgesC1.intersection(edgesC2)
    plausible = True
    H = nx.DiGraph()
    H.add_edges_from(intersectingEdges, data=True)
    if len(H)==0:
        return False
    intersectingNodes = nodesC1.intersection(nodesC2)
    for n in intersectingNodes:
        H.add_node(n)
    for c in nx.weakly_connected_components(H):
        intersectingGComponent = H.subgraph(c).copy()
        if len(intersectingGComponent.edges())%2==0:
            return False
        else:
            for n in intersectingGComponent.nodes():
                if len(intersectingGComponent.in_edges(n))==0:
                    if n in R == True:
                        return False
                if len(intersectingGComponent.out_edges(n))==0:
                    if n in X==True:
                        return False
    return plausible
#############################
#############################


def generateEdgeCycleDict(queue):
    edgeCycleDict = {}
    cycleIDDict = {}
    cycleIDEdgeDict = {}
    cycleIDNodeDict = {}
    visitedEdges = set()
    for i in range(len(queue)):
        c = queue[i]
        cycleIDDict[i] = c
        edgeSet = set()
        nodeSet = set()
        for j in range(len(c)):
            nodeSet.add(c[j])
            if j == len(c)-1:
                e = (c[j], c[0])
            else:
                e = (c[j], c[j+1])
            edgeSet.add(e)
            if e in edgeCycleDict.keys():
                edgeCycleDict[e].add(i)
            else:
                edgeCycleDict[e] = {i}
        cycleIDEdgeDict[i] = edgeSet
        cycleIDNodeDict[i] = nodeSet
        visitedEdges.add(frozenset(edgeSet))
    return edgeCycleDict, cycleIDDict, cycleIDEdgeDict, cycleIDNodeDict, visitedEdges
#############################
############################# 


def getIntersectingCycles(cKey:int, cycleIDEdgeDict:int, edgeCycleDict:dict):
    edgesC = cycleIDEdgeDict[cKey]
    intersectingCycles = set()
    for e in edgesC:
        intersectingCycles = intersectingCycles.union(edgeCycleDict[e])
    return intersectingCycles, edgesC
#############################
#############################


###=============================================================###
###                         Main All                            ###
###=============================================================###

randomGraphs = True
ears = False
fluffles=False
bound=None

if "fluffles" in sys.argv:
    fluffles=True
    randomGraphs=False

if "ears" in sys.argv:
    ears=True
    randomGraphs=False


print("randomGraphs:", randomGraphs, "fluffles:", fluffles, "ears:", ears)
input()

###=====================================================================###
###                    1. Random graphs                                 ###
###=====================================================================###

if randomGraphs==True:
    
    maxSize = 5                                     # (Maximum size is maxSize * l)
    l = int(sys.argv[1])
    a = float(sys.argv[2])
    b = float(sys.argv[3])

    direc= "./Benchmarking/"+str(a)+"to"+str(b)+"DegreeSorted"
    if not os.path.exists(direc):
        os.makedirs(direc)

    bound=None
    timeDict = {}
    cycleDict = {}
    sizeDict = {}
    myList = range(maxSize+1)

    # First: Likelihood
    maxP=7 

    minNoNodes = int((a*1*l)+(b*1*l))
    maxNoNodes = int((a*(len(myList)-1)*l)+(b*(len(myList)-1)*l))
    for k in tqdm(range(1, maxP+1), leave=False, total=maxP):
        p = k/200
        keys = set()
        # Second: Size
        
        for i in tqdm(range(1,len(myList)), leave = False, total=len(myList)-1):
            # Determine sizes
            noMetabolites=int(a*i*l)
            noReactions=int(b*i*l)
            noVertices = noMetabolites + noReactions
            
            if p<=0.015:
                noGraphs = 100000
            elif p>0.015 and p<=0.022:
                noGraphs = 10000
            elif p>0.022 and p <0.03:
                noGraphs = 1000
            else:
                noGraphs = 100
            for j in tqdm(range(noGraphs), leave = False, total=noGraphs):            
                G=nx.algorithms.bipartite.random_graph(n=noMetabolites, m=noReactions, p=p, directed=True)
                X = {n for n,d in G.nodes(data=True) if d['bipartite'] == 0}
                R = set(G) - X
                # assign labels
                
                # Set
                timeStamp = time.time()
                setCounter = 0
                for c in findAllMRChordlessCyclesSet(G, R, X, bound):
                    setCounter += 1
                setTime = time.time()-timeStamp
                
                # List
                timeStamp = time.time()
                listCounter = 0
                for c in findAllMRChordlessCyclesList(G, R, X, bound):
                    listCounter += 1
                listTime = time.time()-timeStamp
                
                # Simple degree
                timeStamp = time.time()
                listDegreeCounter = 0
                for c in findAllMRChordlessCyclesListSimpleDegree(G, R, X, bound):
                    listDegreeCounter+=1
                listDegreeTime = time.time()-timeStamp

                # Simple in-degree
                timeStamp = time.time()
                listInDegreeCounter = 0
                for c in findAllMRChordlessCyclesListInDegree(G, R, X, bound):
                    listInDegreeCounter+=1
                listInDegreeTime = time.time()-timeStamp

                # Simple out-degree
                timeStamp = time.time()
                listOutDegreeCounter = 0
                for c in findAllMRChordlessCyclesListOutDegree(G, R, X, bound):
                    listOutDegreeCounter+=1
                listOutDegreeTime = time.time()-timeStamp

                # R-InDegree
                timeStamp = time.time()
                listRInDegreeCounter = 0
                for c in findAllMRChordlessCyclesListReacInDegree(G, R, X, bound):
                    listRInDegreeCounter += 1
                listRInDegreeTime = time.time()-timeStamp

                # R-OutDegree
                timeStamp = time.time()
                listROutDegreeCounter = 0
                for c in findAllMRChordlessCyclesListReacOutDegree(G, R, X, bound):
                    listROutDegreeCounter += 1
                listROutDegreeTime = time.time()-timeStamp

                # R-SimpleDegree
                timeStamp = time.time()
                listRSymmDegreeCounter = 0
                for c in findAllMRChordlessCyclesListReacSymmDegree(G, R, X, bound):
                    listRSymmDegreeCounter += 1
                listRSymmDegreeTime = time.time()-timeStamp
                
                timeStamp = time.time()
                johnsonCounter = 0
                totalJohnsonCounter = 0
                for c in nx.simple_cycles(G, length_bound=bound):
                    totalJohnsonCounter+=1
                    if checkCycle(c, G, X, R):
                        johnsonCounter+=1
                johnsonTime = time.time()-timeStamp
                
                checkIfCountersMatch(setCounter, listCounter, listDegreeCounter, listROutDegreeCounter, listRSymmDegreeCounter, johnsonCounter)

                assignToTimeDict(timeDict, setCounter, setTime, listTime, listDegreeTime, listInDegreeTime, listOutDegreeTime, listRInDegreeTime, listROutDegreeTime, listRSymmDegreeTime, johnsonTime)
                
                # Assign values to other dictionaries
                if setCounter not in sizeDict.keys():
                    sizeDict[setCounter] = {"setTime": [setTime], "listTime": [listTime], "listDegreeTime": [listDegreeTime], "listInDegreeTime": [listInDegreeTime], "listOutDegreeTime": [listOutDegreeTime],"listRInDegreeTime": [listRInDegreeTime], "listROutDegreeTime": [listROutDegreeTime], "listRSymmDegreeTime":[listRSymmDegreeTime], "johnsonTime": [johnsonTime], "johnsonCounter": [johnsonCounter]}
                else:
                    sizeDict[setCounter]["setTime"].append(setTime)
                    sizeDict[setCounter]["listTime"].append(listTime)
                    sizeDict[setCounter]["listDegreeTime"].append(listDegreeTime)
                    sizeDict[setCounter]["listInDegreeTime"].append(listInDegreeTime)
                    sizeDict[setCounter]["listOutDegreeTime"].append(listOutDegreeTime)                    
                    sizeDict[setCounter]["listRInDegreeTime"].append(listRInDegreeTime)
                    sizeDict[setCounter]["listROutDegreeTime"].append(listROutDegreeTime)
                    sizeDict[setCounter]["listRSymmDegreeTime"].append(listRSymmDegreeTime)
                    sizeDict[setCounter]["johnsonTime"].append(johnsonTime)
                    sizeDict[setCounter]["johnsonCounter"].append(johnsonCounter)
                keys.add(setCounter)
                cycleDict[setCounter]=cycleDict.setdefault(setCounter, [])+ [totalJohnsonCounter]              
                if p>0.02:
                    path = direc+ "/Benchmarking" +str(p) + "min_" + str(minNoNodes) + "max_" +str(maxNoNodes) +".pkl"
                    plotResultsLine(minNoNodes, maxNoNodes, p, direc, timeDict)
                    with open(path, "wb") as file:
                        pickle.dump((timeDict, cycleDict, sizeDict), file)
        path = direc+ "/Benchmarking" +str(p) + "min_" + str(minNoNodes) + "max_" +str(maxNoNodes) +".pkl"
        plotResultsLine(minNoNodes, maxNoNodes, p, direc, timeDict)
        with open(path, "wb") as file:
            pickle.dump((timeDict, cycleDict, sizeDict), file)

###=====================================================================###
###                         Fluffles                                    ###
###=====================================================================###

if fluffles==True:
    maxSize = 5                                     # (Maximum size is maxSize * l)

    l = int(sys.argv[1])
    a = float(sys.argv[2])
    b = float(sys.argv[3])

    direc= "./Benchmarking/"+str(a)+"to"+str(b)+"Fluffles"
    if not os.path.exists(direc):
        os.makedirs(direc)

    bound=None
    timeDict = {}
    cycleDict = {}
    fluffleDict = {}
    overAllTimeDict = {}
    myList = range(maxSize+1)

    # First: Likelihood
    maxP=7 

    for k in tqdm(range(1, maxP+1), leave=False, total=maxP):
        probabilityDict = {}
        p = k/200
        keys = set()
        # Second: Size
        for i in tqdm(range(1,len(myList)), leave = False, total=len(myList)-1):
            sizeDict = {}
            
            # Determine sizes
            noMetabolites=int(a*i*l)
            noReactions=int(b*i*l)
            noVertices = noMetabolites + noReactions
            if i ==1:
                minNoNodes = noMetabolites+noReactions
            elif i==len(myList)-1:
                maxNoNodes = noMetabolites+noReactions
            
            if p<=0.015:
                noGraphs = 100000
            elif p>0.015 and p<=0.022:
                noGraphs = 10000
            elif p>0.022 and p <0.03:
                noGraphs = 1000
            else:
                noGraphs = 100
            for j in tqdm(range(noGraphs), leave = False, total=noGraphs):            
                G=nx.algorithms.bipartite.random_graph(n=noMetabolites, m=noReactions, p=p, directed=True)
                X = {n for n,d in G.nodes(data=True) if d['bipartite'] == 0}
                R = set(G) - X
                # assign labels
                
                # timeStamp = time.time()
                # setCounter = sum(1 for c in findAllMRChordlessCyclesSet(G, R, X, bound))
                # setTime = time.time()-timeStamp
                
                timeStamp = time.time()
                queue = []
                listCounter = 0
                for c in findAllMRChordlessCyclesList(G, R, X, bound):
                    listCounter +=1
                    queue.append(c)
                listTime = time.time()-timeStamp
                
                # # Simple - degree
                # timeStamp = time.time()
                # listDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListSimpleDegree(G, R, X, bound))
                # listDegreeTime = time.time()-timeStamp

                # # R-OutDegree
                # timeStamp = time.time()
                # listROutDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListReacOutDegree(G, R, X, bound))
                # listROutDegreeTime = time.time()-timeStamp

                # # R-SimpleDegree
                # timeStamp = time.time()
                # listRSymmDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListReacSymmDegree(G, R, X, bound))
                # listRSymmDegreeTime = time.time()-timeStamp

                johnsonCounter = 0
                totalJohnsonCounter = 0
                timeStamp = time.time()
                for c in nx.simple_cycles(G, length_bound=bound):
                    totalJohnsonCounter+=1
                    if checkCycle(c, G, X, R):
                        johnsonCounter+=1
                johnsonTime = time.time()-timeStamp
                if johnsonCounter!=listCounter:
                    sys.exit("Counters do not match")
                
                if listCounter>10000:
                    continue
                
                timeStamp=time.time()
                fluffleCounter = listCounter + assembleFluffles(queue=queue, G=G, X=X, R=R)
                fluffleTime = time.time()-timeStamp

                assignToTimeDictEears(timeDict, listCounter, listTime, johnsonTime, fluffleTime)
                
                # Assign values to other dictionaries
                if listCounter not in sizeDict.keys():
                    sizeDict[listCounter] = {"listTime": [listTime], "johnsonTime": [johnsonTime], "johnsonCounter": [johnsonCounter], "fluffleTime": [fluffleTime], "fluffleCounter": [fluffleCounter]}
                else:
                    sizeDict[listCounter]["listTime"].append(listTime)
                    sizeDict[listCounter]["johnsonTime"].append(johnsonTime)
                    sizeDict[listCounter]["johnsonCounter"].append(totalJohnsonCounter)
                    sizeDict[listCounter]["fluffleTime"].append(fluffleTime)
                    sizeDict[listCounter]["fluffleCounter"].append(fluffleCounter)
                keys.add(listCounter)
                cycleDict[listCounter]=cycleDict.setdefault(listCounter, [])+ [totalJohnsonCounter]
                fluffleDict[listCounter]=cycleDict.setdefault(listCounter, [])+ [fluffleCounter]
                probabilityDict[noVertices] = sizeDict

                if p>0.02 and i>0:
                    path = direc+ "/Benchmarking" +str(p) + "min_" + str(minNoNodes) + "max_" +str(maxNoNodes) +".pkl"
                    plotResultsLine(minNoNodes, maxNoNodes, p, direc, timeDict)
                    with open(path, "wb") as file:
                        pickle.dump((timeDict, cycleDict, fluffleDict, overAllTimeDict), file)
        overAllTimeDict[p]=probabilityDict
        path = direc+ "/Benchmarking" +str(p) + "min_" + str(minNoNodes) + "max_" +str(maxNoNodes) +".pkl"
        plotResultsLine(minNoNodes, maxNoNodes, p, direc, timeDict)
        with open(path, "wb") as file:
            pickle.dump((timeDict, cycleDict, fluffleDict, overAllTimeDict), file)


###======================================================================###
###                          Fluffles with MR-ears                       ###
###======================================================================###

# Random MR-EarGraphs
if ears==True:
    maxk = int(sys.argv[1])                 # Max number of ears to be added
    maxl = int(sys.argv[2])                 # Max average length of an ear

    timeDict = {}
    cycleDict = {}
    sizeDict = {}
    fluffleDict = {}
    fluffleMRDict={}
    fluffleJohnsonDict={}
    keys=set()
    direc= "./Benchmarking/Max"+str(maxk)+"EarsMaxLength"+str(maxl)+"RandomMREar"
    if not os.path.exists(direc):
        os.makedirs(direc)

    for k in tqdm(range(1,maxk+1), leave=False, total=maxk+1-1):
        for length in tqdm(range(2,maxl+1),leave=False, total=maxl+1-2):
            if k<=25:
                r=100
            else:
                r=10
            for i in tqdm(range(r), leave=False, total=r):
                G, X, R=randomMREarGraph.mainRandomMREarGraph(k, length)
                # X = {n for n,d in G.nodes(data=True) if d['M'] == True}
                # R = set(G)-X
                if not nx.is_bipartite(G):                    
                    print(G.nodes())
                    print(G.edges())
                    input()
                
                timeStamp = time.time()
                setCounter = sum(1 for c in findAllMRChordlessCyclesSet(G, R, X, bound))
                setTime = time.time()-timeStamp
                
                timeStamp = time.time()
                listCounter = sum(1 for c in findAllMRChordlessCyclesList(G, R, X, bound))
                listTime = time.time()-timeStamp
                
                # Simple - degree
                timeStamp = time.time()
                listDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListSimpleDegree(G, R, X, bound))
                listDegreeTime = time.time()-timeStamp

                # Simple in-degree
                timeStamp = time.time()
                listInDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListInDegree(G, R, X, bound))
                listInDegreeTime = time.time()-timeStamp
                
                # Simple in-degree
                timeStamp = time.time()
                listOutDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListOutDegree(G, R, X, bound))
                listOutDegreeTime = time.time()-timeStamp

                # R-InDegree
                timeStamp = time.time()
                listRInDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListReacInDegree(G, R, X, bound))
                listRInDegreeTime = time.time()-timeStamp

                # R-OutDegree
                timeStamp = time.time()
                listROutDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListReacOutDegree(G, R, X, bound))
                listROutDegreeTime = time.time()-timeStamp

                # R-SimpleDegree
                timeStamp = time.time()
                listRSymmDegreeCounter = sum(1 for c in findAllMRChordlessCyclesListReacSymmDegree(G, R, X, bound))
                listRSymmDegreeTime = time.time()-timeStamp

                # R-SimpleDegree            
                timeStamp = time.time()
                totalJohnsonCounter = sum(1 for c in nx.simple_cycles(G, length_bound=bound) if checkCycle(c,G,X,R))
                johnsonTime = time.time()-timeStamp

                timeStamp = time.time()
                johnsonWOCheckCounter = sum(1 for c in nx.simple_cycles(G, length_bound=bound))
                johnsonWOCheckTime = time.time()-timeStamp

                if johnsonWOCheckCounter!=listCounter:
                    print()
                    print("G is bipartite", nx.is_bipartite(G))
                    print("List counter is", listCounter, "JohnsonCounter is", totalJohnsonCounter, "Total JohnsonCounter is", johnsonWOCheckCounter)
                    sys.exit("Counters do not match!")

                # timeStamp=time.time()
                # fluffleCounter = listCounter + assembleFluffles(queue, G, X, R)
                # fluffleTime = time.time()-timeStamp
                
                checkIfCountersMatch(setCounter, listCounter, listDegreeCounter, listROutDegreeCounter, listRSymmDegreeCounter, totalJohnsonCounter)

                assignToTimeDict(timeDict, setCounter, setTime, listTime, listDegreeTime, listInDegreeTime, listOutDegreeTime, listRInDegreeTime, listROutDegreeTime, listRSymmDegreeTime, johnsonTime)

                if "JohnsonWOCheck" not in timeDict:
                    timeDict["JohnsonWOCheck"]={}
                if setCounter in timeDict["JohnsonWOCheck"].keys():
                    timeDict["JohnsonWOCheck"][setCounter].append(johnsonWOCheckTime)
                else:
                    timeDict["JohnsonWOCheck"][setCounter]=[johnsonWOCheckTime]

                # Assign values to other dictionaries
                if setCounter not in sizeDict.keys():
                    sizeDict[setCounter] = {"setTime": [setTime], "listTime": [listTime], "listDegreeTime": [listDegreeTime], "listInDegreeTime": [listInDegreeTime], "listOutDegreeTime": [listOutDegreeTime],"listRInDegreeTime": [listRInDegreeTime], "listROutDegreeTime": [listROutDegreeTime], "listRSymmDegreeTime":[listRSymmDegreeTime], "johnsonTime": [johnsonTime], "johnsonCounter": [totalJohnsonCounter],
                    "JohnsonWOCheckTime": [johnsonWOCheckTime], "JohnsonWOCheckCounter": [johnsonWOCheckCounter]}
                else:
                    sizeDict[setCounter]["setTime"].append(setTime)
                    sizeDict[setCounter]["listTime"].append(listTime)
                    sizeDict[setCounter]["listDegreeTime"].append(listDegreeTime)
                    sizeDict[setCounter]["listInDegreeTime"].append(listInDegreeTime)
                    sizeDict[setCounter]["listOutDegreeTime"].append(listOutDegreeTime)
                    sizeDict[setCounter]["listRInDegreeTime"].append(listRInDegreeTime)                    
                    sizeDict[setCounter]["listROutDegreeTime"].append(listROutDegreeTime)
                    sizeDict[setCounter]["listRSymmDegreeTime"].append(listRSymmDegreeTime)
                    sizeDict[setCounter]["johnsonTime"].append(johnsonTime)
                    sizeDict[setCounter]["johnsonCounter"].append(totalJohnsonCounter)
                    sizeDict[setCounter]["JohnsonWOCheckTime"].append(johnsonWOCheckTime)
                    sizeDict[setCounter]["JohnsonWOCheckCounter"].append(johnsonWOCheckCounter)
                keys.add(setCounter)
                cycleDict[setCounter]=cycleDict.setdefault(setCounter, [])+ [totalJohnsonCounter]                

                noVertices=len(X)+len(R)
                
                path = direc+ "/Benchmarking_NoEars_" + str(k) + "_MaxEarLength_" + str(length) + ".pkl"                
                with open(path, "wb") as file:
                    pickle.dump((timeDict, cycleDict, sizeDict), file)
                plotResultsLineEars(k,length,direc,timeDict)
        