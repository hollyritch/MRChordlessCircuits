import networkx as nx
from collections import defaultdict
from itertools import product


def chordless_cycle_search_Species(F, B, path, length_bound, G):
    fwBlocked = defaultdict(int)
    bwBlocked = defaultdict(int)
    target = path[0]
    for i in range(1, len(path)):
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
                        for r in Bx:
                            bwBlocked[r] += 1
                    else:
                        for m in Fx:
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


def chordless_cycle_search_Reaction(F:dict, B:dict, path:list, length_bound:int, G:nx.DiGraph, U:nx.DiGraph, overlappingNodes:set):
    ## Reaction is the first node
    
    fwBlocked = defaultdict(int)                                    # will contain only reactions
    bwBlocked = defaultdict(int)                                    # will contain only metabolites‚
    target = path[0]                                                # Target is a reaction
    current = path[1]
    for i in range(len(path)):                                      # Initializing
        a = path[i]
        if i==2:                                                    # Reaction case 
            Ba=U.predecessors(a)
            for m in Ba:
                if m in overlappingNodes:
                    m_in = str(m)+"_in"
                    m_out = str(m)+"_out"
                    bwBlocked[m_in] +=1
                    bwBlocked[m_out] +=1
                else:
                    if m == current:
                        Ba_m = B[a]
                        if m in Ba_m:
                            bwBlocked[m] +=1                                    # Block m from being visited (only unvisited metabolites can be visited)
                    else:
                        bwBlocked[m] +=1
        else:   
            Fa = U.successors(a)                                              # Species case
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
                if sum(1 for y in Fx if y in path)==1: 
                    yield path + [x]                    # yield the path and stop.....you can' go any further, otherwise there would be non MR-chordfree cycles
            else:                                   
                if G.nodes[x]["Type"]=="Species":   # Block forward in the case of a metabolite
                    if type(x)==str:
                        x_= int(x.split("_")[0])
                    else:
                        x_ = x
                    Fx_ = U.successors(x_)            # We cannot have x=m so we can do the lookup in the U.succ
                    for y in Fx_:
                        fwBlocked[y] += 1
                    # for y in Fx:
                    #     fwBlocked[y] += 1
                else:                               # Again for reactions we need to make sure that they also 
                                                    # update the overlapping metabolites, but be carefull with the current....for this one only the updates according to the current oriented network           
                    Bx=U.predecessors(x)
                    for m in Bx:
                        if m in overlappingNodes:
                            m_in = str(m)+"_in"
                            m_out = str(m)+"_out"
                            bwBlocked[m_in] +=1
                            bwBlocked[m_out] +=1
                        else:
                            if m == current:
                                Bx_m = B[x]
                                if m in Bx_m:           # here we make sure that there is really an edge in the current 
                                                        # out-directed network
                                    bwBlocked[m] +=1    # if this is not the case....we don't want the update (we take 
                                                        # care of these cycles in the other out-direcred network)
                            else:
                                bwBlocked[m] +=1        
                    # Bx = B[x]
                    # for y in Bx:                      # Block backward in case of a reaction
                    #     bwBlocked[y] += 1
                path.append(x)                          # Either way, append x to the path add all forward oriented vertices 
                                                        # to the stack
                stack.append(iter(Fx))
                break
        else:                                           # Take off 
            stack.pop() 
            z = path.pop()
            if G.nodes[z]["Type"]=="Species":
                if type(z)==str:
                    z_=int(z.split("_")[0])
                else:
                    z_ = z
                Fz_=U.successors(z_)
                for y in Fz_:
                    fwBlocked[y] -= 1
            else:                                       # The reaction do not need to be modified 
                Bz = U.predecessors(z)
                for y in Bz:
                    if y in overlappingNodes:
                        y_in = str(y)+"_in"
                        y_out = str(y)+"_out"
                        bwBlocked[y_in] -= 1
                        bwBlocked[y_out] -= 1
                    else:
                        if y == current:
                            Bz_y = B[z]
                            if y in Bz_y:           # here we make sure that there is really an edge in the current 
                                                    # out-directed network
                                bwBlocked[y] -=1    # if this is not the case....we don't want the update (we take 
                                                    # care of these cycles in the other out-direcred network)
                        else:
                            bwBlocked[y] -=1        
#############################
#############################


def chordless_cycle_search(F:dict, B:dict, path:list, length_bound, G:nx.DiGraph, U:nx.DiGraph, overlappingMetabolites:set):
    if G.nodes[path[0]]["Type"]=="Species":
        yield from chordless_cycle_search_Species(F, B, path, length_bound, G)
    else:
        yield from chordless_cycle_search_Reaction(F, B, path, length_bound, G, U, overlappingMetabolites)
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
        Bc = B.subgraph(c)
        Fcc = Bcc = None
        for S in stems(Fc, v):
            if Fcc is None:
                Fcc = nx.algorithms.cycles._NeighborhoodCache(Fc)
                Bcc = nx.algorithms.cycles._NeighborhoodCache(Bc)
            yield from chordless_cycle_search(Fcc, Bcc, S, bound, F, F, set())
        components.extend(c for c in nx.strongly_connected_components(F.subgraph(c - {v})) if len(c) > 3)
#############################
#############################


def findNodeSpecificMRChordlessCycles(F, U, x, bound, overlappingMetabolites:set):
    B = F.reverse(copy=True)
    
    # Directed stems look like (u -> v -> w), so we use the product of
    # predecessors of v with successors of v.
    # Findet praktisch Triangles oder alle möglichen u->v->w Pfade von allen Vorgängern zu allen nach Folgern von v
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
                yield from chordless_cycle_search(Fcc, Bcc, S, bound, F, U, overlappingMetabolites)
            break
#############################
#############################


G = nx.DiGraph()

G.add_node("M1")
G.add_node("M2")
G.add_node("M3")
G.add_node("M4")
# G.add_node("M5")
# G.add_node("M6")
# G.add_node("M7")
# G.add_node("M8")
# G.add_node("M9")
# G.add_node("M10")

G.add_node("R1")
G.add_node("R2")
G.add_node("R3")
G.add_node("R4")
# G.add_node("R5")
# G.add_node("R6")
# G.add_node("R7")
# G.add_node("R8")
# G.add_node("R9")
# G.add_node("R10")

G.nodes["M1"]["Type"] = "Species"
G.nodes["R1"]["Type"] = "Reaction"
G.nodes["M2"]["Type"] = "Species"
G.nodes["R2"]["Type"] = "Reaction"
G.nodes["M3"]["Type"] = "Species"
G.nodes["R3"]["Type"] = "Reaction"
G.nodes["M4"]["Type"] = "Species"
G.nodes["R4"]["Type"] = "Reaction"


G.add_edge("M1", "R1")
G.add_edge("R1", "M2")
G.add_edge("R1", "M3")

G.add_edge("M2", "R2")
G.add_edge("M1", "R2")
G.add_edge("R2", "M3")

G.add_edge("M3", "R3")
G.add_edge("R3", "M1")



for c in findNodeSpecificMRChordlessCycles(G, "M3", None, G, set()):
    print(c)