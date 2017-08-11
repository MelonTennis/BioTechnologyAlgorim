import MedianString2_4
import MotifEnum2_1
import random

# String spelled by genome path problem
def genomePath(pattern):
    k = len(pattern[0])
    res = pattern[0]
    for i in range(1, len(pattern)):
        res = res + pattern[i][-1]
    return res

# solve the overlap graph problem
def overlapGrapth(Dna):
    res = []
    for i in range(0, len(Dna)):
        for j in range(0, len(Dna)):
            if i != j and Dna[i][1:] == Dna[j][:-1]:
                res.append(Dna[i] + ' -> ' + Dna[j])
    return res

# construct deBruijin graph
def deBruijin(k, Dna):
    diction = dict()
    for i in range(0, len(Dna) - k):
        key = Dna[i: i+k]
        if key not in diction:
            diction[key] = []
        diction[key].append(Dna[i+1: i+1+k])
    res = []
    for m in diction:
        cur = m + " -> "
        for i in range(0, len(diction[m])):
            cur = cur + diction[m][i]
            if i != len(diction[m]) - 1:
                cur = cur + ", "
        res.append(cur)
    return sorted(res)

# construct k-mer deBruijin graph
def kmerDebrujin(patterns):
    diction = dict()
    for pattern in patterns:
        key = pattern[:-1]
        if key not in diction:
            diction[key] = []
        diction[key].append(pattern[1:])
    res = []
    for m in diction:
        cur = m + " -> "
        for i in range(0, len(diction[m])):
            cur = cur + diction[m][i]
            if i != len(diction[m]) - 1:
                cur = cur + ", "
        res.append(cur)
    return sorted(res)


# LinkedList Node for form cycle
class Node:
    def __init__(self, str):
        self.val = str
        self.next = None
        self.visit = False
        self.adj = []

class absNode:
    def __init__(self, v):
        self.val = v
        self.next = None

# form a cycle, return the entrance of the cycle, as a node
def formCycle(edges, start, available, nodes):
    #print edges[start.val]
    nextNode = Node(edges[start.val][0])
    if nextNode.val in nodes:
        nextNode = nodes[nextNode.val]
    else:
        nodes[nextNode.val] = nextNode
    firstPass = True
    #if start.next != None:
    #    start.adj.append(start.next)
    cur = start
    while nextNode.val != start.val or firstPass:
        firstPass = False
        if cur.next == None:
            cur.next = nextNode
            follow = cur.next
        else:
            cur.adj.append(nextNode)
            follow = cur.adj[-1]
        edges[cur.val].remove(nextNode.val)
        if len(edges[cur.val]) == 0:
            edges.pop(cur.val, None)
            if cur in available:
                available.remove(cur)
        else:
            if cur not in available:
                available.append(cur)
        cur = follow
        nextNode = Node(edges[cur.val][0])
        if nextNode.val in nodes:
            nextNode = nodes[nextNode.val]
        else:
            nodes[nextNode.val] = nextNode
        if nextNode.val == start.val:
            edges[cur.val].remove(nextNode.val)
            if len(edges[cur.val]) > 0 and cur not in available:
                available.append(cur)
            elif len(edges[cur.val]) == 0:
                edges.pop(cur.val, None)
                if cur in available:
                    available.remove(cur)
        #print edges
        #print "available: ", [x.val for x in available]
    if cur.next == None:
        cur.next = start
    else:
        if start not in cur.adj:
            cur.adj.append(start)
    if len(available) != 0:
        return available[0]
    else:
        return start

# combine two cycles
def combine(cycle, cycle1):
    cycle1.adj.append(cycle1.next)
    while cycle.val != cycle1.val:
        if len(cycle.adj) != 0:
            return combine(cycle.adj[0], cycle1)
        else:
            cycle = cycle.next
    cycle1.next = cycle.next
    return cycle1

# solve the Eulerian Cycle Problem
def EeulerianCycle(graph):
    edges = dict()
    available = []
    nodes = dict()
    absNodes = dict()
    for line in graph:
        start = line.strip().split(' ')[0]
        if start not in edges:
            edges[start] = []
        adj = line.strip().split(' ')[2].split(',')
        for nei in adj:
            edges[start].append(nei)
    #print edges
    begin = Node(edges.keys()[0])
    nodes[begin.val] = begin
    cycle = formCycle(edges, begin, available, nodes)
    while len(available) != 0:
        #print "available: ", [x.val for x in available]
        #print "edges: ", edges
        assert (len(edges) != 0 and len(available) != 0)
        cycle = formCycle(edges, available[0], available, nodes)
    #print "nodes", nodes.keys()
    #print [(x.val, x.next.val, [m.val for m in x.adj]) for x in [nodes[y] for y in nodes]]
    assert(len(available) == 0)
    cur = cycle
    available.append(cur)
    head = None
    while len(available) != 0:
        current = available[0]
        cur = current
        firstPass = True
        if cur.val in absNodes:
            head = absNodes[cur.val][0]
        else:
            head = absNode(cur.val)
            absNodes[cur.val] = [head]
        record = head.next
        while cur.val != current.val or firstPass:
            firstPass = False
            temp = cur.next
            if temp == None :
                print cur.val
                print "edges: ", edges
                print "available: ", [x.val for x in available]
                print "should not get temp == None"
            if len(cur.adj) != 0:
                follow = cur.adj.pop(0)
                cur.next = follow
            else:
                cur.next = None
            if cur.next != None:
                if cur not in available:
                    available.append(cur)
            else:
                if cur in available:
                    available.remove(cur)
            head.next = absNode(temp.val)
            if temp.val in absNodes:
                absNodes[temp.val].append(head.next)
            else:
                absNodes[temp.val] = [head.next]
            head = head.next
            cur = temp
        #head.next = absNode(current.val)
        head.next = record

    res = ""
    head = absNodes[cycle.val][0]
    firstPass = True
    while head.next != None:
        if firstPass:
            res = res + head.val
        else:
            res = res + "->" + head.val
        head = head.next
        firstPass = False
    res = res + "->" + cycle.val
    return res
'''

class Node:
    def __init__(self, val):
        self.val = val
        self.next = None

def EeulerianCycle(graph):
    edges = dict()
    available = dict()
    nodes = dict()
    for line in graph:
        start = line.strip().split(' ')[0]
        if start not in edges:
            edges[start] = []
        adj = line.strip().split(' ')[2].split(',')
        for nei in adj:
            edges[start].append(nei)
    current = edges.keys()[0]
    firstKey = 0
    for key in edges.keys():
        if len(edges[key]) > 1:
            current = key
            firstKey = key
            print key
            break
    #print edges
    #head = Node(current)
    #nodes[current] = [head]
    available[current] = edges[current]
    while len(available) != 0:
        current = available.keys()[0]
        firstPass = True
        if current in nodes:
            head = nodes[current][0]
        else:
            head = Node(current)
            nodes[current] = [head]
        record = head.next
        while head.val != current or firstPass:
            firstPass = False
            nextVal = edges[head.val].pop(0)
            nextNode = Node(nextVal)
            if nextVal in nodes:
                nodes[nextVal].append(nextNode)
            else:
                nodes[nextVal] = [nextNode]
            head.next = nextNode
            if len(edges[head.val]) == 0:
                edges.pop(head.val, None)
                if head.val in available:
                    available.pop(head.val, None)
            else:
                if head.val not in available:
                    available[head.val] = edges[head.val]
            head = head.next
        #head.next = Node(current)
        head.next = record
    print "nodes", [(x, len(nodes[x])) for x in nodes.keys()]

    res = ""
    head = nodes[nodes.keys()[0]][0]
    firstPass = True
    while head.next != None:
        if firstPass:
            res = res + head.val
        else:
            res = res + "->" + head.val
        head = head.next
        firstPass = False
    res = res + "->" + firstKey
    return res
'''

# form a path in graph
def formPath(edges, start, available, nodes, first, last):
    #print edges[start.val]
    nextNode = Node(edges[start.val][0])
    if nextNode.val in nodes:
        nextNode = nodes[nextNode.val]
    else:
        nodes[nextNode.val] = nextNode
    firstPass = True
    #if start.next != None:
    #    start.adj.append(start.next)
    cur = start
    while nextNode.val != start.val or firstPass:
        firstPass = False
        if cur.next == None:
            cur.next = nextNode
            follow = cur.next
        else:
            cur.adj.append(nextNode)
            follow = cur.adj[-1]
        edges[cur.val].remove(nextNode.val)
        if len(edges[cur.val]) == 0:
            edges.pop(cur.val, None)
            if cur in available:
                available.remove(cur)
        else:
            if cur not in available:
                available.append(cur)
        cur = follow
        if cur.val == last and (cur.val not in edges or len(edges[cur.val]) == 0):
            break
        nextNode = Node(edges[cur.val][0])
        if nextNode.val in nodes:
            nextNode = nodes[nextNode.val]
        else:
            nodes[nextNode.val] = nextNode
        if nextNode.val == start.val:
            edges[cur.val].remove(nextNode.val)
            if len(edges[cur.val]) > 0 and cur not in available:
                available.append(cur)
            elif len(edges[cur.val]) == 0:
                edges.pop(cur.val, None)
                if cur in available:
                    available.remove(cur)
        #print edges
        #print "available: ", [x.val for x in available]
    if cur.next == None:
        cur.next = start
    else:
        if start not in cur.adj:
            cur.adj.append(start)
    return start

# return a Eulerian Path in graph
def EeulerianPath(graph):
    edges = dict()
    available = []
    nodes = dict()
    absNodes = dict()
    pathHead = None
    pathEnd = None
    indegree = dict()
    outdegree = dict()
    for line in graph:
        start = line.strip().split(' ')[0]
        if start not in edges:
            edges[start] = []
        adj = line.strip().split(' ')[2].split(',')
        if start not in outdegree:
            outdegree[start] = len(adj)
        for nei in adj:
            if nei not in indegree:
                indegree[nei] = 1
            else:
                indegree[nei] = indegree[nei] + 1
            edges[start].append(nei)
    remain = [x for x in indegree.keys() if x not in outdegree or indegree[x] > outdegree[x]]
    #print remain
    if len(remain) != 0:
        pathEnd = remain.pop()
    assert(len(remain) == 0)
    odd = [x for x in outdegree.keys() if x not in indegree.keys()]
    if len(odd) != 0:
        pathHead = odd.pop()
    else:
        odd = [x for x in outdegree.keys() if indegree[x] < outdegree[x]]
        pathHead = odd.pop(0)
        if len(odd) != 0:
            pathEnd = odd.pop()
    print pathHead, pathEnd
    #print outdegree[pathHead], indegree[pathEnd], outdegree[pathEnd]
    begin = Node(pathHead)
    nodes[begin.val] = begin
    cycle = formPath(edges, begin, available, nodes, pathHead, pathEnd)
    while len(available) != 0:
        #print "available: ", [x.val for x in available]
        #print "edges: ", edges
        assert (len(edges) != 0 and len(available) != 0)
        cycle1 = formPath(edges, available[0], available, nodes, pathHead, pathEnd)
    #print "nodes", nodes.keys()
    #print [(x.val, x.next.val, [m.val for m in x.adj]) for x in [nodes[y] for y in nodes]]
    assert(len(available) == 0)
    cur = cycle
    available.append(cur)
    head = None
    while len(available) != 0:
        current = available[0]
        cur = current
        firstPass = True
        if cur.val in absNodes:
            head = absNodes[cur.val][0]
        else:
            head = absNode(cur.val)
            absNodes[cur.val] = [head]
        record = head.next
        while cur.val != current.val or firstPass:
            firstPass = False
            temp = cur.next
            if temp == None :
                print cur.val
                print "edges: ", edges
                print "available: ", [x.val for x in available]
                print "should not get temp == None"
            if len(cur.adj) != 0:
                follow = cur.adj.pop(0)
                cur.next = follow
            else:
                cur.next = None
            if cur.next != None:
                if cur not in available:
                    available.append(cur)
            else:
                if cur in available:
                    available.remove(cur)
            head.next = absNode(temp.val)
            if temp.val in absNodes:
                absNodes[temp.val].append(head.next)
            else:
                absNodes[temp.val] = [head.next]
            head = head.next
            cur = temp
        #head.next = absNode(current.val)
        head.next = record
    res = ""
    head = absNodes[pathHead][0]
    firstPass = True
    while head.next != None:
        if firstPass:
            res = res + head.val
        else:
            res = res + "->" + head.val
        head = head.next
        firstPass = False
    return res

# construct dna from patterns
def kPattern(k, patterns):
    edges = dict()
    available = []
    nodes = dict()
    absNodes = dict()
    pathHead = None
    pathEnd = None
    indegree = dict()
    outdegree = dict()
    for line in patterns:
        start = line.strip()[0:k-1]
        end = line.strip()[1:]
        if start not in edges:
            edges[start] = [end]
        else:
            edges[start].append(end)
        if start not in outdegree:
            outdegree[start] = 1
        else:
            outdegree[start] = outdegree[start] + 1
        if end not in indegree:
            indegree[end] = 1
        else:
            indegree[end] = indegree[end] + 1
    remain = [x for x in indegree.keys() if x not in outdegree or indegree[x] > outdegree[x]]
    # print remain
    if len(remain) != 0:
        pathEnd = remain.pop()
    assert (len(remain) == 0)
    odd = [x for x in outdegree.keys() if x not in indegree.keys()]
    if len(odd) != 0:
        pathHead = odd.pop()
    else:
        odd = [x for x in outdegree.keys() if indegree[x] < outdegree[x]]
        pathHead = odd.pop(0)
        if len(odd) != 0:
            pathEnd = odd.pop()
    print pathHead, pathEnd
    print edges
    # print outdegree[pathHead], indegree[pathEnd], outdegree[pathEnd]
    begin = Node(pathHead)
    nodes[begin.val] = begin
    cycle = formPath(edges, begin, available, nodes, pathHead, pathEnd)
    while len(available) != 0:
        # print "available: ", [x.val for x in available]
        # print "edges: ", edges
        assert (len(edges) != 0 and len(available) != 0)
        cycle1 = formPath(edges, available[0], available, nodes, pathHead, pathEnd)
    # print "nodes", nodes.keys()
    # print [(x.val, x.next.val, [m.val for m in x.adj]) for x in [nodes[y] for y in nodes]]
    assert (len(available) == 0)
    cur = cycle
    available.append(cur)
    head = None
    while len(available) != 0:
        current = available[0]
        cur = current
        firstPass = True
        if cur.val in absNodes:
            head = absNodes[cur.val][0]
        else:
            head = absNode(cur.val)
            absNodes[cur.val] = [head]
        record = head.next
        while cur.val != current.val or firstPass:
            firstPass = False
            temp = cur.next
            if temp == None:
                print cur.val
                print "edges: ", edges
                print "available: ", [x.val for x in available]
                print "should not get temp == None"
            if len(cur.adj) != 0:
                follow = cur.adj.pop(0)
                cur.next = follow
            else:
                cur.next = None
            if cur.next != None:
                if cur not in available:
                    available.append(cur)
            else:
                if cur in available:
                    available.remove(cur)
            head.next = absNode(temp.val)
            if temp.val in absNodes:
                absNodes[temp.val].append(head.next)
            else:
                absNodes[temp.val] = [head.next]
            head = head.next
            cur = temp
        # head.next = absNode(current.val)
        head.next = record
    res = ""
    head = absNodes[pathHead][0]
    firstPass = True
    while head.next != None:
        if firstPass:
            res = res + head.val
        else:
            res = res + head.val[-1]
        head = head.next
        firstPass = False
    return res

# Solve the String Reconstruction from Read-Pairs Problem.
def readPairs(k, d, patterns):
    edges = dict()
    available = []
    nodes = dict()
    absNodes = dict()
    pathHead = None
    pathEnd = None
    indegree = dict()
    outdegree = dict()
    for line in patterns:
        front = line.split('|')[0].strip()
        back = line.split('|')[1].strip()
        start = (front[:-1], back[:-1])
        end = (front[1:], back[1:])
        if start not in edges:
            edges[start] = [end]
        else:
            edges[start].append(end)
        if start not in outdegree:
            outdegree[start] = 1
        else:
            outdegree[start] = outdegree[start] + 1
        if end not in indegree:
            indegree[end] = 1
        else:
            indegree[end] = indegree[end] + 1
    remain = [x for x in indegree.keys() if x not in outdegree or indegree[x] > outdegree[x]]
    # print remain
    if len(remain) != 0:
        pathEnd = remain.pop()
    assert (len(remain) == 0)
    odd = [x for x in outdegree.keys() if x not in indegree.keys()]
    if len(odd) != 0:
        pathHead = odd.pop()
    else:
        odd = [x for x in outdegree.keys() if indegree[x] < outdegree[x]]
        pathHead = odd.pop(0)
        if len(odd) != 0:
            pathEnd = odd.pop()
    print pathHead, pathEnd
    #print edges
    # print outdegree[pathHead], indegree[pathEnd], outdegree[pathEnd]
    begin = Node(pathHead)
    nodes[begin.val] = begin
    cycle = formPath(edges, begin, available, nodes, pathHead, pathEnd)
    while len(available) != 0:
        # print "available: ", [x.val for x in available]
        # print "edges: ", edges
        assert (len(edges) != 0 and len(available) != 0)
        cycle1 = formPath(edges, available[0], available, nodes, pathHead, pathEnd)
    # print "nodes", nodes.keys()
    # print [(x.val, x.next.val, [m.val for m in x.adj]) for x in [nodes[y] for y in nodes]]
    assert (len(available) == 0)
    cur = cycle
    available.append(cur)
    head = None
    while len(available) != 0:
        current = available[0]
        cur = current
        firstPass = True
        if cur.val in absNodes:
            head = absNodes[cur.val][0]
        else:
            head = absNode(cur.val)
            absNodes[cur.val] = [head]
        record = head.next
        while cur.val != current.val or firstPass:
            firstPass = False
            temp = cur.next
            if temp == None:
                print cur.val
                print "edges: ", edges
                print "available: ", [x.val for x in available]
                print "should not get temp == None"
            if len(cur.adj) != 0:
                follow = cur.adj.pop(0)
                cur.next = follow
            else:
                cur.next = None
            if cur.next != None:
                if cur not in available:
                    available.append(cur)
            else:
                if cur in available:
                    available.remove(cur)
            head.next = absNode(temp.val)
            if temp.val in absNodes:
                absNodes[temp.val].append(head.next)
            else:
                absNodes[temp.val] = [head.next]
            head = head.next
            cur = temp
        # head.next = absNode(current.val)
        head.next = record
    res = ""
    head = absNodes[pathHead][0]
    firstPass = True
    count = 0
    tail = ""
    while head.next != None:
        if firstPass:
            res = res + head.val[0]
            tail = tail + head.val[1]
        else:
            res = res + head.val[0][-1]
            tail = tail + head.val[1][-1]
        head = head.next
        firstPass = False
    ans = res[0:k+d]+tail
    return ans

# Solve the Contig Generation Problem.
def contig(patterns):
    nodes = dict()
    indegree = dict()
    outdegree = dict()
    edges = dict()
    for line in patterns:
        start = line[:-1]
        end = line[1:]
        if start not in edges:
            edges[start] = [end]
        else:
            edges[start].append(end)
        if start not in outdegree:
            outdegree[start] = 1
        else:
            outdegree[start] = outdegree[start] + 1
        if end not in indegree:
            indegree[end] = 1
        else:
            indegree[end] = indegree[end] + 1
    inonly = [x for x in indegree.keys() if x not in outdegree.keys()]
    outonly = [x for x in outdegree.keys() if x not in indegree.keys()]
    for m in inonly:
        outdegree[m] = 0
    for m in outonly:
        indegree[m] = 0
    assert(len(indegree.keys()) == len(outdegree.keys()))
    paths = []
    for key in indegree.keys():
        if not indegree[key] == outdegree[key] == 1:
            if outdegree[key] > 0:
                for adj in edges[key]:
                    NoBpath = [key, adj]
                    while indegree[adj] == outdegree[adj] == 1:
                        adj = edges[adj][0]
                        NoBpath.append(adj)
                    paths.append(NoBpath)
    #print paths
    res = []
    res = [makeRes(x) for x in paths]
    return res

# trun to right form
def makeRes(path):
    res = ""
    first = True
    for p in path:
        if first:
            res = res + p
        else:
            res = res + p[-1]
        first = False
    return res




infile = '/Users/yijia/Documents/current work/Bio/dataset_38427_5.txt'
with open(infile) as f:
    Dna = f.readlines()
Dna = [x.strip() for x in Dna[0:]]
res = contig(Dna)
#print res
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38427_5_res.txt', 'w')
for p in res:
    outfile.write(p)
    outfile.write('\n')
outfile.close()