import sys

def topological(adjList):
    list1 = []
    indegree = dict()
    outdegree = dict()
    candidate = set()
    for key in adjList.keys():
        if key not in outdegree:
            outdegree[key] = len(adjList[key])
        else:
            outdegree[key] = outdegree[key] + len(adjList[key])
        if key not in indegree:
            indegree[key] = 0
        for x in adjList[key]:
            if x in indegree:
                indegree[x] = indegree[x] + 1
            else:
                indegree[x] = 1
            if x not in outdegree:
                outdegree[x] = 0
    #print "in", indegree
    #print "out", outdegree
    for m in indegree.keys():
        if indegree[m] == 0:
            candidate.add(m)
    while len(candidate) != 0:
        a = candidate.pop()
        list1.append(a)
        if a not in adjList.keys():
            continue
        for i in range(len(adjList[a])-1, -1, -1):
            b = adjList[a][i]
            indegree[b] = indegree[b] - 1
            if indegree[b] == 0:
                candidate.add(b)
            adjList[a].remove(b)
            if len(adjList[a]) == 0:
                adjList.pop(a, None)
    if len(adjList) != 0:
        return "not DAG"
    else:
        return list1

def longestPath(adj, source, sink, weight):
    nodes = set()
    for k in adj.keys():
        nodes.add(k)
        for bb in adj[k]:
            nodes.add(bb)
    s = dict()
    path = dict()
    for n in nodes:
        s[n] = -10000000000
        path[n] = []
    s[source] = 0
    topo = topological(adj)
    breakpoint = 0
    for i in range(0, len(topo)):
        if topo[i] == source:
            breakpoint = i
            break
    topo = topo[breakpoint:]
    ancess = set()
    while len(topo) != 0:
        b = topo.pop(0)
        for j in ancess:
            str = j + '->' + b
            if str in weight:
                s[b] = max(s[b], weight[str] + s[j])
                if s[b] == weight[str] + s[j]:
                    path[b].append(j)
        ancess.add(b)
    longpath = [sink]
    cur = sink
    while cur != source:
        longpath.append(path[cur][-1])
        cur = path[cur][-1]
    res = longpath[-1]
    longpath.reverse()
    for p in longpath[1:]:
        res = res + '->' + p
    print s[sink]
    return res

def DPChange(money, coins):
    miniCoin = [0 for x in range(0, money+1)]
    for m in range(1, money+1):
        miniCoin[m] = sys.maxint
        for i in range(0, len(coins)):
            if m >= coins[i]:
                if miniCoin[m - coins[i]] + 1 < miniCoin[m]:
                    miniCoin[m] = miniCoin[m - coins[i]] + 1
    return miniCoin[money]

def ManhattanTourist(n, m, down, right):
    #print down
    #print right
    s = [[0 for i in range(0, m+1)] for x in range(0, n+1)]
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
        #print s[i]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + right[0][j-1]
        #print s[0][j]
    #print s
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])
        #print s[i]
    return s[n][m]

def LCSBackTrack(v, w):
    back = [["" for i in range(0, len(w))] for i in range(0, len(v))]
    s = [[0 for i in range(0, len(w)+1)] for i in range(0, len(v)+1)]
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            s[i][j] = max(s[i-1][j], s[i][j-1])
            if v[i-1] == w[j-1]:
                s[i][j] = max(s[i-1][j], s[i-1][j-1] + 1, s[i][j-1])
            if s[i][j] == s[i-1][j]:
                back[i-1][j-1] = "d"
            elif s[i][j] == s[i][j-1]:
                back[i-1][j-1] = "r"
            elif s[i][j] == s[i-1][j-1] + 1 and v[i-1] == w[j-1]:
                back[i-1][j-1] = "w"
    return back

def OutputLCS(backtrack, v, ii ,jj):
    res = ""
    i = ii
    j = jj
    while i >= 0 and j >= 0:
        if backtrack[i][j] == "d":
             i = i-1
        elif backtrack[i][j] == "r":
            j = j-1
        else:
            res = v[i] + res
            i = i-1
            j = j-1
    return res



infile = '/Users/yijia/Documents/current work/Bio/dataset_38467_7.txt'
with open(infile) as f:
    input = f.readlines()
adj = dict()
weight = dict()
start = input[0].strip()
end = input[1].strip()
for p in input[2:]:
    arr = p.strip().split('->')
    arrw = p.strip().split(':')
    key = arr[0].strip()
    value = arr[1].strip().split(":")
    if key not in adj:
        adj[key] = [value[0]]
    else:
        adj[key] = adj[key] + [value[0]]
    weight[arrw[0]] = int(arrw[1])

#print weight
res = longestPath(adj, start, end, weight)
#print longPath

#n = int(input[0].strip().split(' ')[0])
#m = int(input[0].strip().split(' ')[1])
#right = []
#down = []
#for i in range(1, 1+n):
#    down.append([int(x) for x in input[i].strip().split(" ")])
#for j in range(2+n, 2*n + 3):
#    right.append([int(x) for x in input[j].strip().split(" ")])

#money = int(input[0].strip())
#coins = [int(x) for x in input[1].strip().split(",")]
#coins = [24,13,12,7,5,3,1]
#money = 8074
#v = input[0].strip()
#w = input[1].strip()
#back = LCSBackTrack(v, w)
#res = OutputLCS(back, v,len(v)-1, len(w)-1)
#print res
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38467_7_res.txt', 'w')
#for i in range(0, len(res)):
#    p = str(res[i])
#    outfile.write(p)
    #if i != len(res) - 1:
    #    outfile.write(', ')
outfile.write(str(res))
outfile.close()