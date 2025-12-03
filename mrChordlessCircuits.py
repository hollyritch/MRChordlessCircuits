from collections import defaultdict
import networkx as nx
from itertools import product

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
                    if x=="M2":
                        print(x, fwBlocked["R2"])
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
                    print("Unfortunately....it does not work")
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


def findAllMRChordlessCycles(F, bound):
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

G=nx.DiGraph()

G.add_node("M1")
G.add_node("M2")
G.add_node("M3")
G.add_node("M4")
G.add_node("M5")

G.add_node("R1")
G.add_node("R2")
G.add_node("R3")
G.add_node("R4")
G.add_node("R5")

G.nodes["M1"]["Type"] = "Species"
G.nodes["R1"]["Type"] = "Reaction"
G.nodes["M2"]["Type"] = "Species"
G.nodes["R2"]["Type"] = "Reaction"
G.nodes["M3"]["Type"] = "Species"
G.nodes["R3"]["Type"] = "Reaction"
G.nodes["M4"]["Type"] = "Species"
G.nodes["R4"]["Type"] = "Reaction"
G.nodes["M5"]["Type"] = "Species"
G.nodes["R5"]["Type"] = "Reaction"

G.add_edge("M1", "R1")
G.add_edge("R1", "M3")

G.add_edge("M2", "R1")
G.add_edge("M2", "R2")
G.add_edge("R2", "M5")

G.add_edge("M3", "R3")
G.add_edge("R3", "M2")

G.add_edge("M4", "R4")
G.add_edge("R4", "M1")

G.add_edge("M5", "R5")
G.add_edge("R5", "M4")


for c in findAllMRChordlessCycles(G, None):
    print(c)


print("Now all cycles")
for c in nx.simple_cycles(G):
    print(c)