import random
import sys
import math
import numpy as np

# 11.6.2
def farestFirstTravel(data, k, d):
    centers = rndChoose(data, 1)
    while len(centers) < k:
        point = MaxPoint(data, centers)
        centers.append(point)
    return centers

def rndChoose(data, k):
    res = []
    for i in range(0, k):
        res.append(data[i])
    return res

def MaxPoint(data, centers):
    point = []
    maxDist = -sys.maxint
    for da in data:
        minDist = sys.maxint
        for center in centers:
            d = math.sqrt(sum([pow(center[i] - da[i], 2) for i in range(0, len(da))]))
            if d < minDist:
                minDist = d
        if minDist > maxDist:
            maxDist = minDist
            point = list(da)
    return point

# 11.7.3
def squareErrorDistortion(data, centers):
    res = 0
    for j in range(0, len(data)):
        dmin = sys.maxint
        for center in centers:
            d = sum([pow(center[i] - data[j][i], 2) for i in range(0, len(data[j]))])
            if d < dmin:
                dmin = d
        res += dmin
    return res/len(data)

def Lloyd(k, dd, data):
    centers = rndChoose(data, k)
    clusters = [[] for i in range(0, k)]
    change = True
    for da1 in data:
        dmin1 = sys.maxint
        curC1 = -1
        for m1 in range(0, len(centers)):
            d1 = math.sqrt(sum([pow(centers[m1][i] - da1[i], 2) for i in range(0, len(da1))]))
            if d1 < dmin1:
                dmin1 = d1
                curC1 = m1
        clusters[curC1].append(da1)
    while change:
        change = False
        for s in range(0, len(clusters)):
            centers[s] = [float(sum(col))/len(col) for col in zip(*clusters[s])]
        for m in range(0, len(clusters)):
            for r in range(len(clusters[m])-1, -1, -1):
                da = clusters[m][r]
                curC = m
                minDist = sys.maxint
                for c in range(0, len(centers)):
                    d = math.sqrt(sum([pow(centers[c][i] - da[i], 2) for i in range(0, len(da))]))
                    if d < minDist:
                        minDist = d
                        curC = c
                if curC != m:
                    change = True
                    clusters[curC].append(da)
                    clusters[m].remove(da)
        centers1= [[] for w in range(0, len(centers))]
        for s in range(0, len(clusters)):
            centers1[s] = [float(sum(col))/len(col) for col in zip(*clusters[s])]
        if sorted(centers1) != sorted(centers):
            change = True
    return centers

# 11.13.7
def EM(k, m, beta, data):
    centers = rndChoose(data, k)
    for iter in range(0, 100):
        centers = HM(data, centers, beta, k)
    return centers

def distance(a, b):
    return math.sqrt(sum([pow(a[ii] - b[ii], 2) for ii in range(0, len(a))]))

def HM(data, centers, beta, k):
    matrix = [[0 for i in range(0, len(data))] for j in range(0, k)]
    for j in range(0, len(data)):
        all = 0
        for i in range(0, k):
            all += math.exp(-beta*distance(centers[i], data[j]))
        for i in range(0, k):
            matrix[i][j] = math.exp(-beta*distance(centers[i], data[j]))/all
    centers = []
    for line in matrix:
        center = [0 for w in range(0, len(data[0]))]
        for i in range(0, len(line)):
            center[i] += sum([line[ii]*data[i][ii] for ii in range(0, len(line))])
        sumCenter = sum(center)
        center = [center[ii]/sumCenter for ii in range(0, len(center))]
        centers.append(center)
    return centers

# 11.14.7
def HierarchicalClustering(n, data):
    clusters = [[i] for i in range(0, n)]
    while len(clusters) > 1:
        ci, cj = twoCloseCluster(clusters, data)
        newC = ci + cj
        clusters.remove(ci)
        clusters.remove(cj)
        clusters.append(newC)
        temp = ""
        for c in newC:
            temp += str(c+1) + " "
        print temp[:-1]

def twoCloseCluster(clusters, data):
    ci = []
    cj = []
    dmin = sys.maxint
    for i in range(0, len(clusters)):
        for j in range(i+1, len(clusters)):
            dvag = 0.0
            for ii in range(len(clusters[i])):
                for jj in range(len(clusters[j])):
                    dvag += data[clusters[i][ii]][clusters[j][jj]]
            dvag = dvag/(len(clusters[i]) * len(clusters[j]))
            if dvag < dmin:
                dmin = dvag
                ci = clusters[i]
                cj = clusters[j]
    return ci, cj

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
k = int(input[0].strip().split(" ")[0])
# d = int(input[0].strip().split(" ")[1])
# beta = float(input[1].strip())
# centers = [[float(x) for x in r.strip().split(" ")] for r in input[1:1+k]]
data = [[float(x) for x in r.strip().split(" ")] for r in input[1:]]
# res = squareErrorDistortion(data, centers)
# res = farestFirstTravel(data, k, d)
# res = Lloyd(k, d, data)
# res = EM(k, d, beta, data)
HierarchicalClustering(k, data)
# print res
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
# for r in res:
#     for rr in range(0, d):
#         outfile.write(str("%.3f"%r[rr]))
#         if rr != d-1:
#             outfile.write(" ")
#     outfile.write("\n")

