import pickle
from matplotlib import pyplot as plt
import numpy as np
import os

def plotTimeResultsLine(path:str, sizeDict:dict, k, l):
    fig, ax = plt.subplots()
    fig.set_figheight(15)
    fig.set_figwidth(15)
    barLabel=["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:cyan", "tab:purple", "tab:olive", "tab:pink"]
    i = 0
    n=32
    totalJohnsonMeans, totalJohnsonStds = [], []
    totalStd=[]
    mrCount = sorted(list(sizeDict.keys()))
    for mr in mrCount:
        totalJohnsonMeans.append(np.mean(list(c for c in sizeDict[mr]["johnsonCounter"])))
        totalJohnsonStds.append(np.std(list(c for c in sizeDict[mr]["johnsonCounter"])))    
    #ax.plot(mrChordlessMeans, ratio, markersize=2,  color="tab:blue")
    ax.errorbar(x=mrCount, y=totalJohnsonMeans, yerr=totalJohnsonStds, elinewidth=0.25, color="tab:green")
    #ax.legend(loc='upper left', ncols=1, fontsize=n)
    ax.set_xlabel("MR-chordless circuits", fontsize=n)
    ax.set_ylabel("Elementary circuits", fontsize=n)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.yticks(fontsize=n)
    plt.xticks(fontsize=n)
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(path + "/MRChordlessVsJohnson_Line_MaxEears"+str(k) + "maxAvEarLength"+str(l)+".png")
    plt.close()


    # # Second plot
    # fig, ax = plt.subplots()
    # fig.set_figheight(15)
    # fig.set_figwidth(15)
    # #ax.plot(mrChordlessMeans, ratio, markersize=2,  color="tab:blue")
    # ax.errorbar(x=fluffleCounts, y=mrChordlessMeans, yerr=mrChordlessStds, fmt="o", markersize=2, capsize=3, elinewidth=0.25, label="Johnson", color="tab:green")
    # #ax.legend(loc='upper left', ncols=1, fontsize=n)
    # ax.set_xlabel("Fluffle count", fontsize=n)
    # ax.set_ylabel("MR-chordless circuits", fontsize=n)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # plt.yticks(fontsize=n)
    # plt.xticks(fontsize=n)
    # if not os.path.exists(path):
    #     os.makedirs(path)
    # plt.savefig(path + "/MRChordlessVsFluffles_Line_MaxEears"+str(k) + "maxAvEarLength"+str(l)+".png")
    return
#############################
#############################

path = "./Benchmarking/Max50EarsMaxLength20RandomMREar/Benchmarking_NoEars_22_MaxEarLength_16.pkl"

k=path.split("NoEars_")[1].split("MaxEarLength")[0]
l=path.split("MaxEarLength_")[1].split(".pkl")[0]
newPath = path.split("Benchmarking_NoEars")[0]
with open(path, "rb") as file:
    data = pickle.load(file)

sizeDict = data[2]

plotTimeResultsLine(newPath, sizeDict, k, l)
