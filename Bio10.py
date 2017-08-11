import sys

# 10.1.7
def neighbourJoin(n, distMap):
    global m
    # print n, distMap
    if n == 2:
        edges = {}
        # print distMap
        for key in distMap.keys():
            edges[key] = dict()
            for key2 in distMap.keys():
                if key2 != key:
                    edges[key][key2] = distMap[key][key2]
        return edges
    totalDist = dict()
    for key in distMap.keys():
        totalDist[key] = sum(distMap[key].values())
    Dstar = neighbourMatrix(distMap)
    ci, cj = findMin(Dstar)
    delta = float(totalDist[ci] - totalDist[cj])/(n-2)
    limbi = float(distMap[ci][cj] + delta)*0.5
    limbj = float(distMap[ci][cj] - delta)*0.5
    if m not in distMap.keys():
        distMap[m] = dict()
    for key in distMap.keys():
        if key != m and key != ci and key != cj:
            distMap[key][m] = float(distMap[key][ci] + distMap[key][cj] - distMap[ci][cj])*0.5
            distMap[m][key] = distMap[key][m]
    distMap.pop(ci, None)
    distMap.pop(cj, None)
    for key in distMap:
        if ci in distMap[key].keys():
            distMap[key].pop(ci, None)
        if cj in distMap[key].keys():
            distMap[key].pop(cj, None)
    # print len(distMap.keys())
    m += 1
    edges = neighbourJoin(n-1, distMap)
    m -= 1
    if m not in edges.keys():
        edges[m] = dict()
    edges[m][ci] = limbi
    if ci not in edges.keys():
        edges[ci] = dict()
    edges[ci][m] = limbi
    edges[m][cj] = limbj
    if cj not in edges.keys():
        edges[cj] = dict()
    edges[cj][m] = limbj

    #print edges
    return edges

def neighbourMatrix(distMap):
    n = len(distMap.keys())
    res = dict()
    totalDist = dict()
    for key in distMap.keys():
        totalDist[key] = sum(distMap[key].values())
    for key in distMap.keys():
        res[key] = dict()
        for key2 in distMap.keys():
            if key2 != key:
                res[key][key2] = (n-2)*distMap[key][key2] - totalDist[key] - totalDist[key2]
    #print res
    return res

def findMin(Dstar):
    ci = -1
    cj = -1
    minValue = sys.maxint
    for i in Dstar.keys():
        for j in Dstar.keys():
            if i != j and Dstar[i][j] < minValue:
                minValue = Dstar[i][j]
                ci = i
                cj = j
    # print ci, cj
    return ci, cj

def matrixToMap(distance):
    distMap = {}
    for i in range(0, len(distance)):
        if i not in distMap.keys():
            distMap[i] = dict()
        for j in range(0, len(distance)):
            distMap[i][j] = distance[i][j]
    #print distMap
    return distMap

def changeForm(map):
    res = []
    for i in map.keys():
        for j in map[i].keys():
            if map[i][j] != 0:
                res.append(str(i) + "->" + str(j) + ":" + "%.3f"%map[i][j])
    return res

# 10.3.10
def smallParsimony(T, characters):
    Tag = {}
    S= {}
    minValue = sys.maxint
    for v in T.keys():
        Tag[v] = 0
        if len(T[v].keys()) == 0:
            Tag[v] = 1
        for k in ['A', 'G', 'C', 'T']:
            if k not in S.keys():
                S[k] = {}
            if characters[v] == k:
                S[k][v] = 0
            else:
                S[k][v] = sys.maxint
    for t in T.keys():
        cnt = 0
        for child in T.keys():
            if Tag[child] == 1:
                cnt == 1
            else:
                break
        if cnt == 2:
            Tag[t] = 1
        for symbol in ['A', 'G', 'C', 'T']:
            phi = 1
            if symbol == t:
                phi = 0
            S[t][v] = min(S[t][list(S.keys())[0]] + phi, S[t][list(S.keys())[1]] + phi)
            minValue = min(minValue, S[t][v])
    return minValue

alpha = ['A', 'C', 'G', 'T']
length = 0
class edge(object):
    def __init__(self):
        self.id = 0
        self.score = 0

class node(object):
    def __init__(self):
        self.edges = []
        self.res = ""
        self.leaf = False
        self.score = []

def solveDP(tree, cur):
    global length
    global alpha
    vertex = tree[cur]
    if vertex.leaf:
        return
    vertex.score = [[0 for i in range(0, length)] for j in range(0, 4)]
    for edge in vertex.edges:
        edge_id = edge.id
        next_vertex = tree[edge_id]
        if len(next_vertex.score) == 0:
            solveDP(tree, edge_id)

        for i in range(0, length):
            for j in range(0, 4):
                curMin = sys.maxint
                if next_vertex.leaf:
                    curMin = 0
                    if j == alpha.index(next_vertex.res[i]):
                        curMin = 1
                else:
                    for k in range(0, 4):
                        add = 0
                        if j != k:
                            add = 1
                        curMin = min(curMin, next_vertex.score[k][i] + add)
                vertex.score[j][i] += curMin

def setTree(tree, cur, prescore):
    global length
    global alpha
    vertex = tree[cur]
    if vertex.leaf == True:
        return
    for edge in vertex.edges:
        next_id = edge.id
        next_vertex = tree[next_id]
        score = 0
        if len(next_vertex.res) == 0:
            for i in range(0, length):
                j = alpha.index(vertex.res[i])
                idx = 0
                add = 1
                if j == 0:
                    add = 0
                minScore = next_vertex.score[0][i] + add
                for k in range(1, 4):
                    add = 1
                    if j == k:
                        add = 0
                    if next_vertex.score[k][i] + add < minScore:
                        minScore = next_vertex.score[k][i] + add
                        idx = k
                next_vertex.res += alpha[idx]
                add = 1
                if j == idx:
                    add = 0
                score += add
            edge.score = score
            setTree(tree, next_id, score)
        elif next_vertex.leaf:
            for i in range(0, length):
                add = 0
                if vertex.res[i] != next_vertex.res[i]:
                    add = 1
                score += add
            edge.score = score
            next_vertex.edges[0].score = score
        else:
            edge.score = score

def smallParsimony(tree, root, score):
    global length
    solveDP(tree, root)
    score = 0
    vertex_root = tree[root]
    for i in range(0, length):
        idx = 0
        minScore = vertex_root.score[0][i]
        for j in range(1, 4):
            if vertex_root.score[j][i] < minScore:
                idx = j
                minScore = vertex_root.score[j][i]
        score += minScore
        vertex_root.res += alpha[idx]
    setTree(tree, root, -1)
    return tree

cur_num = 0
def buildTree(n, adj):
    global length
    global cur_num
    tree = {}
    for couple in adj:
        begin_id = int(couple[0])
        end_id = -1
        if begin_id not in tree.keys():
            tree[begin_id] = node()
        if couple[1][0] <= '9' and couple[1][0] >= '0':
            end_id = int(couple[1])
            if end_id not in tree.keys():
                tree[end_id] = node()
        else:
            if length == 0:
                length = len(couple[1])
            end_id = cur_num
            cur_num += 1
            tree[end_id] = node()
            cur_leaf = tree[end_id]
            cur_leaf.leaf = True
            cur_leaf.res = couple[1]
        edge_begin = edge()
        edge_end = edge()
        edge_begin.id = end_id
        edge_end.id = begin_id
        tree[begin_id].edges.append(edge_begin)
        tree[end_id].edges.append(edge_end)
    return tree

# 10.4.6
def nearestNeighbour(a, b, adj):
    res1 = list(adj)
    res2 = list(adj)
    candidate1 = set()
    candidate2 = set()
    for edge in adj:
        if (edge[0] == a and edge[1] != b) or (edge[0] == b and edge[1] != a):
            if edge[0] == a:
                candidate1.add(edge[1])
            else:
                candidate2.add(edge[1])
            res1.remove(edge)
            res2.remove(edge)
        if (edge[1] == a and edge[0] != b) or (edge[1] == b and edge[0] != a):
            if edge[1] == a:
                candidate1.add(edge[0])
            else:
                candidate2.add(edge[0])
            res1.remove(edge)
            res2.remove(edge)
    candidate1 = list(candidate1)
    candidate2 = list(candidate2)
    res1.append((candidate1[0], a))
    res1.append((a, candidate1[0]))
    res2.append((candidate1[1], b))
    res2.append((b, candidate1[1]))
    res1.append((candidate2[1], b))
    res1.append((b, candidate2[1]))
    res2.append((candidate2[1], a))
    res2.append((a, candidate2[1]))
    res1.append((candidate2[0], a))
    res1.append((a, candidate2[0]))
    res2.append((candidate2[0], b))
    res2.append((b, candidate2[0]))
    res1.append((candidate1[1], b))
    res1.append((b, candidate1[1]))
    res2.append((candidate1[0], a))
    res2.append((a, candidate1[0]))
    return res1, res2

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
n = int(input[0].strip())
# k = int(input[1].strip())
# pa = [int(x) for x in input[0].strip().split(" ")]
adj = [(m[0], m[1]) for m in [x.strip().split('->') for x in input[1:]]]
#lines = [[int(r) for r in line.strip().split(" ")] for line in input[1:]]
#m = n
# lines = matrixToMap(lines)
# res1, res2 = nearestNeighbour(pa[0], pa[1], adj)
tree = buildTree(n, adj)
for r in tree.keys():
    if len(tree[r].edges) == 2:
        root = r
        break
res = smallParsimony(tree, root, 0)
for key in res.keys():
    e = res[key]
    print key, e.res, e.score
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
# for r in res:
#     outfile.write(str(r[0]))
#     outfile.write('->')
#     outfile.write(str(r[1]))
#     outfile.write("\n")
# outfile.write("\n")
# for r in res2:
#     outfile.write(str(r[0]))
#     outfile.write('->')
#     outfile.write(str(r[1]))
#     outfile.write("\n")

