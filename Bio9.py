import sys

# 9.2.12
def readLine(lines):
    adj = dict()
    for line in lines:
        adjs = line.split('->')
        head = int(adjs[0])
        end = int(adjs[1].split(':')[0])
        score = int(adjs[1].split(':')[1])
        if head not in adj.keys():
            adj[head] = dict()
        if end not in adj[head].keys():
            adj[head][end] = score
        if end not in adj.keys():
            adj[end] = dict()
        if head not in adj[end].keys():
            adj[end][head] = score
    return adj

def constructMatrix(k, lines):
    adj = readLine(lines)
    matrix = findPath(adj)
    # for line in matrix:
    #     print line
    res = [matrix[i][0:k] for i in range(0, k)]
    return res

def findPath(adj):
    size = len(adj.keys())
    res = [[-1 for ii in range(0, size)] for jj in range(0, size)]
    for i in range(0, size):
        res[i][i] = 0
        getScore(i, i, adj, res, 0, size)
    print res
    return res

def getScore(i, next, adj, res, pre, size):
    for j in adj[next].keys():
        if res[i][j] == -1:
            res[i][j] = adj[next][j] + pre
            res[j][i] = res[i][j]
            for p in adj[j].keys():
                if res[i][p] == -1:
                    getScore(i, j, adj, res, res[i][j], size)

# 9.3.11
def LimbLength(n, j, distance):
    res = sys.maxint
    for i in range(0, n):
        for k in range(0, n):
            if i != j and k != j:
                value = (distance[i][j] + distance[j][k] - distance[i][k])/2.0
                res = min(value, res)
    return int(res)

# 9.4.6
N = 29
def additivePhylogeny(distance, n):
    global N
    if n == 2:
        paths = dict()
        paths[0] = {1: distance[0][1]}
        paths[1] = {0: distance[1][0]}
        return paths
    limb = LimbLength(n, n-1, distance[0:n][0:n])
    for i in range(0, n-1):
        distance[i][n-1] -= limb
        distance[n-1][i] -= limb
    for k in range(1, n-1):
        if distance[0][k] == distance[0][n-1] + distance[n-1][k]:
            break
    x = distance[0][n-1]
    paths = additivePhylogeny(distance, n-1)
    i = 0
    path = searchPath(i, k, paths, [])
    cur = 0
    left = x
    for idx in range(0, len(path)-1):
        p = path[idx]
        cur += p[1]
        if cur > x:
            paths[i][N] = left
            if N not in paths.keys():
                paths[N] = dict()
            paths[N][n-1] = limb
            paths[N][i] = left
            if n-1 not in paths.keys():
                paths[n-1] = dict()
            paths[n-1][N] = limb
            paths[N][path[idx+1][0]] = p[1] - left
            if path[idx+1][0] not in paths.keys():
                paths[path[idx+1][0]] = dict()
            paths[path[idx + 1][0]] = p[1] - left
            N += 1
            break
        else:
            i = p[0]
            left -= p[1]
    return paths

def searchPath(begin, end, paths, visit):
    res = []
    for p in paths.keys():
        if p == begin:
            if end in paths[p].keys():
                visit.append(end)
                res.append([end, paths[begin][end]])
                return res
            else:
                for q in paths[p].keys():
                    if q not in visit:
                        visit.append(q)
                        res = res + searchPath(q, end, paths, visit)
            break
    return res

# 9.6.8
def UPGMA(distance, n):
    N = n
    cluster = {}
    paths = {}
    age = {}
    distmap = {}
    for i in range(0, len(distance)):
        distmap[i] = {}
        for j in range(0, len(distance)):
            if i != j:
                distmap[i][j] = distance[i][j]
    for i in range(0, n):
        age[i] = 0.0
        cluster[i] = [i]
    while len(cluster) > 1:
        ci, cj = closeCluster(distmap)
        # print ci, cj
        cnew = cluster[ci] + cluster[cj]
        cluster[N] = cnew
        if ci not in paths.keys():
            paths[ci] = {}
        if cj not in paths.keys():
            paths[cj] = {}
        if N not in paths.keys():
            paths[N] = {}
        age[N] = float(distmap[ci][cj]) / 2.0
        paths[ci][N] = age[N] - age[ci]
        paths[cj][N] = age[N] - age[cj]
        paths[N][ci] = age[N] - age[ci]
        paths[N][cj] = age[N] - age[cj]
        distmap = editDistance(ci, cj, N, distmap, cluster)
        cluster.pop(ci, None)
        cluster.pop(cj, None)
        N += 1
    root = list(cluster.keys())[0]
    # for v in paths.keys():
    #     for w in paths.keys():
    #         if v != w:
    #             paths[v][w] = abs(age[v] - age[w])
    #             paths[w][v] = paths[v][w]
    return paths

def closeCluster(cluster, distmap):
    dist = float(sys.maxint)
    ci = -1
    cj = -1
    for s in cluster.keys():
        for e in cluster.keys():
            if s != e:
                d = 0.0
                for s1 in cluster[s]:
                    for e1 in cluster[e]:
                        if s1 != e1 and s1 in distmap.keys():
                            d += distmap[s1][e1]
                d = float(d)/(len(cluster[s])*len(cluster[e]))
                #print d
                if d < dist:
                    dist = d
                    ci = s
                    cj = e
    return ci, cj

def closeCluster(distmap):
    dist = float(sys.maxint)
    ci = -1
    cj = -1
    for i in distmap.keys():
        for j in distmap.keys():
            if i != j and distmap[i][j] < dist:
                dist = distmap[i][j]
                ci = i
                cj = j
    return ci, cj

def editDistance(ci, cj, N, distmap, cluster):
    distmap[N] = {}
    for cm in distmap.keys():
        if cm != ci and cm != cj and cm != N:
            distmap[N][cm] = (distmap[ci][cm] * len(cluster[ci]) + distmap[cj][cm] * len(cluster[cj]))/float(len(cluster[N]))
            distmap[cm][N] = distmap[N][cm]
    distmap.pop(ci, None)
    for d in distmap.keys():
        distmap[d].pop(ci, None)
    distmap.pop(cj, None)
    for d in distmap.keys():
        distmap[d].pop(cj, None)
    return distmap

def changeForm(map):
    res = []
    for i in map.keys():
        for j in map[i].keys():
            if map[i][j] != 0:
                res.append(str(i) + "->" + str(j) + ":" + "%.3f"%map[i][j])
    return res


infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
n = int(input[0].strip())
# k = int(input[1].strip())
lines = [[int(r) for r in line.strip().split(" ")] for line in input[1:]]
res = UPGMA(lines, n)
res = changeForm(res)
print res
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
for r in res:
    outfile.write(r)
    outfile.write("\n")



