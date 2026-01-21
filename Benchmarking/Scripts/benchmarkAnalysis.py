from matplotlib import pyplot as plt
import pickle
import numpy as np
from scipy.special import expit
import os 


def gompertz(x, a, b, c):
    """
    Gompertz function:
    - a: asymptote (max or starting value)
    - b: horizontal shift
    - c: growth rate (positive for increasing, negative for decreasing)
    """
    return a * np.exp(-b * np.exp(-c * x))
#############################
#############################

def plotTimeResultsLineEars(path:str, timeDict:dict, pngPath:str, j:int):
    fig, ax = plt.subplots()
    fig.set_figheight(15)
    fig.set_figwidth(15)
    barLabel=["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink", "tab:brown", "tab:gray"]
    i = 0
    n=32
    for algo, sizeDict in timeDict.items():
        col = barLabel[i]
        if algo == "MR-chordlessSet":
            name="MR set"
            i+=1
            if j!=1:
                continue
        if algo == "MR-chordlessList":
            name="MR list"
            i+=1        
            if j!=2:
                continue
        if algo=="MR-ChordlessListSimpleDegree":
            name="MR list D"        
            i+=1         
            if j!=3:
                continue
        elif algo=="listInDegree":
            name="MR list ID"
            i+=1
            if j!=4:
                continue
        elif algo=="listOutDegree":
            name="MR list OD"
            i+=1
            if j!=5:
                continue
        elif algo=="MR-ChordlessListRInDegreeTime":
            name="MR list RID"
            i+=1
            if j!=7:
                continue
        elif algo=="MR-ChordlessListROutDegreeTime":
            name="MR list ROD"            
            i+=1
            if j!=8:
                continue
        elif algo=="MR-ChordlessListRSymmDegreeTime":
            name="MR list RD"            
            i+=1
            if j!=6:
                continue
        elif algo=="Johnson":
            name="Johnson+"            
            i+=1
        elif algo=="JohnsonWOCheck":
            name="Johnson"            
            i+=1
        keys = sorted(list(sizeDict.keys()))
        values = []
        
        for l in keys:
            values.append(np.mean(sizeDict[l]))
        ax.plot(keys, values, linewidth=1, label=name, color=col)

    ax.legend(loc='upper left', ncols=1, fontsize=n)
    ax.set_xlabel("Count", fontsize=n)
    ax.set_ylabel("Time in seconds", fontsize=n)
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(pngPath)
    plt.close()
    return
#############################
#############################


def plotTimeResultsLine(path:str, timeDict:dict, pngPath:str):
    fig, ax = plt.subplots()
    fig.set_figheight(15)
    fig.set_figwidth(15)
    barLabel=["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink", "tab:brown", "tab:gray"]
    i = 0
    n=32
    for algo, sizeDict in timeDict.items():
        col = barLabel[i]
        if algo == "MR-chordlessSet":
            name="MR set"
            i+=1
        if algo == "MR-chordlessList":
            name="MR list"
            i+=1        
        if algo=="MR-ChordlessListSimpleDegree":
            name="MR list D"        
            i+=1         
        elif algo=="MR-ChordlessListInDegree":
            name="MR list ID"
            i+=1            
        elif algo=="MR-ChordlessListOutDegree":
            name="MR list OD"
            i+=1
        elif algo=="MR-ChordlessListRInDegreeTime":
            name="MR list RID"
            i+=1
        elif algo=="MR-ChordlessListROutDegreeTime":
            name="MR list ROD"            
            i+=1
        elif algo=="MR-ChordlessListRSymmDegreeTime":
            name="MR list RD"            
            i+=1
        elif algo=="Johnson":
            name="Johnson+"            
            i+=1
        elif algo=="JohnsonWOCheck":
            name="Johnson"            
            i+=1
        keys = sorted(list(sizeDict.keys()))
        values = []
        
        for l in keys:
            values.append(np.mean(sizeDict[l]))
        ax.plot(keys, values, linewidth=1, label=name, color=col)

    ax.legend(loc='upper left', ncols=1, fontsize=n)
    ax.set_xlabel("Count", fontsize=n)
    ax.set_ylabel("Time in seconds", fontsize=n)
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(pngPath)
    plt.close()
    return
#############################
#############################

def plotResultsLine(cycleDict:dict, path:str):
    identifier = path.split("/")[2]
    newPath = "./Benchmarking/"+path.split("/")[2]+"/" +identifier+ "_sizes.png"
    fig, ax = plt.subplots()
    fig.set_figheight(15)
    fig.set_figwidth(15)
    print(path)
    print(newPath)
    barLabel=["tab:blue", "orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink"]
    i = 0
    n=32
    x = list(sorted(cycleDict.keys()))[1:]
    y_mean = []
    y_std = []
    for l in x:
        y_mean.append(np.mean(list(l/v for v in cycleDict[l] if v!=0)))
        y_std.append(np.std(list(l/v for v in cycleDict[l] if v!=0)))
    
    #y_smooth = savgol_filter(y_mean, window_length=5, polyorder=2)
    
    ax.errorbar(x, y_mean, yerr=y_std, fmt="o", markersize=2, capsize=3, elinewidth=0.25)
    ax.set_xscale("log")
    ax.set_xlabel("Count of MR-chordless circuits", fontsize = n)
    ax.set_ylabel("Ratio MR-chordless/elementary circuits", fontsize = n)
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    plt.savefig(newPath)
    plt.close()
    return
#############################
#############################


def sigmoid(x, L, k, x0):
    return L - L / (1 + expit(-k * (x - x0)))

# pathList = ["1.0to2.0/Benchmarking0.025min_30max_150.pkl", "1.0to2.0DegreeSorted/Benchmarking0.025min_30max_150.pkl", "1.0to1.0DegreeSorted/Benchmarking0.025min_30max_150.pkl", "1.0to1.0/Benchmarking0.025min_30max_150.pkl", "2.0to1.0DegreeSorted/Benchmarking0.025min_30max_150.pkl", "2.0to1.0/Benchmarking0.025min_30max_150.pkl",
#             "1.0to3.0/Benchmarking0.025min_32max_160.pkl", "3.0to1.0/Benchmarking0.025min_32max_160.pkl", "1.0to4.0/Benchmarking0.025min_30max_150.pkl",
#             "4.0to1.0/Benchmarking0.035min_30max_150.pkl", "1.0to5.0/Benchmarking0.035min_30max_150.pkl", "5.0to1.0/Benchmarking0.035min_30max_150.pkl", "1.0to1.0Fluffles/Benchmarking0.015min_30max_150.pkl", "1.0to2.0Fluffles/Benchmarking0.015min_30max_150.pkl", "2.0to1.0Fluffles/Benchmarking0.015min_30max_150.pkl",
#             "Max50EarsMaxLength20RandomMREar/Benchmarking_NoEars_11_MaxEarLength_3.pkl"]
pathList = ["1.0to2.0DegreeSorted/Benchmarking0.025min_30.0max_150.0.pkl", "1.0to1.0DegreeSorted/Benchmarking0.025min_30.0max_150.0.pkl", "2.0to1.0DegreeSorted/Benchmarking0.025min_30.0max_150.0.pkl", "Max50EarsMaxLength20RandomMREar/Benchmarking_NoEars_25_MaxEarLength_20.pkl"]

for p in pathList:
    print(p)
    path ="./Benchmarking/"+p
    with open(path, "rb") as file:
        data = pickle.load(file)

    cycleDict = data[1]    

    timeDict = data[0]
    overallTimeDict = data[2]
    identifier = p.split("/")[0]
    newPath = "./Benchmarking/"+identifier + "/newPlots"
    if "NoEars" in p:
        noEars = int(p.split("NoEars_")[1].split("_MaxEarLength")[0])
        maxLength = int(p.split("MaxEarLength_")[1].split(".")[0])
        pngPath=newPath + "/MRChordlessVsJohnson_Line_MaxNoEars_"+str(noEars)+"_MaxLength_"+str(maxLength)
        for j in range(1,9):
            newPNGPath=newPath+"/EarsTime"+str(j)
            plotTimeResultsLineEars(newPath, timeDict, newPNGPath, j)            
    else:
        minNoNodes = int(p.split("max_")[0].split("_")[1].split(".")[0])
        maxNoNodes = int(p.split("max_")[1].split(".")[0].split(".")[0])
        prob = float(p.split("Benchmarking")[1].split("min")[0])
        pngPath= newPath + "/MRChordlessVsJohnson_Line_"+str(prob)+"min_" + str(minNoNodes) + "max_" +str(maxNoNodes)+".png"
        plotTimeResultsLine(newPath, timeDict, pngPath)
    plotResultsLine(cycleDict, path)
    
    
