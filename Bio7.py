#7.4.4
def GreedySorting(nums, outfile):
    dist = 0
    for i in range(0, len(nums)):
        if abs(nums[i]) != i+1:
            nums = kSort(nums, i)
            temp = [str(n) if n < 0 else "+" + str(n) for n in nums]
            temp = str(temp)[1:-1].replace(",", "")
            temp = str(temp)[1:-1].replace("'", "")
            outfile.write("("+temp+")")
            outfile.write('\n')
            dist = dist + 1
        if nums[i] == -(i+1):
            nums = kSort(nums, i)
            temp = [str(n) if n < 0 else "+" + str(n) for n in nums]
            temp = str(temp)[1:-1].replace(",", "")
            temp = str(temp)[1:-1].replace("'", "")
            outfile.write("("+temp+")")
            outfile.write('\n')
            dist = dist + 1
    return dist

def kSort(nums, k):
    if abs(nums[k]) != k+1:
        record = k
        while(record < len(nums) and abs(nums[record]) != k+1):
            record = record + 1
        temp = [-nums[m] for m in range(k, record+1)]
        temp.reverse()
        nums = nums[:k] + temp + nums[record+1:]
    else:
        nums[k] = abs(nums[k])
    #print nums
    return nums

def calculateBreakPoint(nums):
    res = 0
    for i in range (1, len(nums)):
        if nums[i] - nums[i-1] != 1:
            res = res + 1
    if nums[0] != 1:
        res = res + 1
    if nums[-1] != len(nums):
        res = res + 1
    return res

# 7.12.4
def chrometocycle(chrome):
    nodes = [0 for i in range(0, len(chrome)*2)]
    for j in range(1, len(chrome)+1):
        i = chrome[j-1]
        if i > 0:
            nodes[2*j-1-1] = 2*i - 1
            nodes[2*j - 1] = 2*i
        else:
            nodes[2*j-1-1] = -2*i
            nodes[2*j-1] = -2*i - 1
    return nodes

# 7.12.5
def cyclechromosome(nodes):
    chromosome = [0 for i in range(0, len(nodes)/2)]
    for j in range(1, len(nodes)/2+1):
        if nodes[2*j-1-1] < nodes[2*j-1]:
            chromosome[j-1] = nodes[2*j-1]/2
        else:
            chromosome[j-1] = -nodes[2*j-1-1]/2
    return chromosome

def listToStr(chrom):
    temp = [str(n) if n < 0 else "+" + str(n) for n in chrom]
    temp = str(temp)[1:-1].replace(",", "")
    temp = str(temp)[1:-1].replace("'", "")
    return "(" + temp + ")"

# 7.12.7
def colorededges(P):
    edges = set()
    for chromo in P:
        nodes = chrometocycle(chromo)
        for j in range(1, len(chromo)+1):
            if 2*j < len(nodes):
                edges.add((nodes[2*j-1], nodes[2*j]))
            else:
                edges.add((nodes[2 * j - 1], nodes[0]))
    return edges

def blackedges(P):
    edges = set()
    p = sum(P, [])
    for i in range(1, len(p)+1):
        edges.add((2*i-1, 2*i))
    return edges

# 7.12.8
def graphtogenome(genomeGraph):
    p = list()
    cycles = constructCycle(genomeGraph)
    for nodes in cycles:
        chrom = cyclechromosome(nodes)
        chrom = [str(n) if n < 0 else "+" + str(n) for n in chrom]
        chrom = str(chrom)[1:-1].replace(",", "")
        chrom = str(chrom)[1:-1].replace("'", "")
        p.append("(" + chrom + ")")
    return p

def constructCycle(graph):
    i = 1
    res = []
    while i < len(graph):
        record = i-1
        while (i < len(graph) and i%2 == 0 and (graph[i-1] == graph[i] + 1 or graph[i-1] == graph[i] - 1))or(i < len(graph) and i%2 == 1):
            i = i + 1
        # begin = graph[i]
        # record = i
        # while i < len(graph) and graph[i] != begin - 1 and graph[i] != begin + 1 or (i - record)%2 != 1:
        #     i = i+1
        res.append(graph[record:i])
        i = i + 1
    res = [m[-1:]+m[0:-1] for m in res]
    #print res
    return res

# 7.9.4
def twoBreakDistance(P, Q):
    block = len(P)
    edges1 = colorededges(P)
    edges2 = colorededges(Q)
    cycle = findCycle(edges1, edges2)
    return block - cycle

def findCycle(edges1, edges2):
    res = 0
    length = len(edges1)
    edges1 = list(edges1)
    edges2 = list(edges2)
    while len(edges1) != 0 and len(edges2) != 0:
        head = edges1[0]
        cur = edges1[0]
        edges1.remove(head)
        while cur[1] != head[0]:
            for i in edges2:
                if i[0] == cur[1] or i[1] == cur[1]:
                    if(i[0] == cur[1]):
                        cur = i
                    else:
                        cur = (i[1], i[0])
                    edges2.remove(i)
                    break
            for j in edges1:
                if j[0] == cur[1] or j[1] == cur[1]:
                    if j[0] == cur[1]:
                        cur = j
                    else:
                        cur = (j[1], j[0])
                    edges1.remove(j)
                    break
        res = res + 1
    #print length - res
    return length - res

# 7.13.2
def twoBreakOnGenomeGraph(graph, i, i1, j, j1):
    if (i, i1) in graph:
        idx = graph.index((i, i1))
        graph[idx] = (i1, j1)
    else:
        idx = graph.index((i1, i))
        graph[idx] = (i1, j1)
    if(j,j1) in graph:
        idx = graph.index((j, j1))
        graph[idx] = (j, i)
    else:
        idx = graph.index((j1, j))
        graph[idx] = (j, i)
    return graph

def twoBreakOnGenome(P, i, i1, j, j1):
    graph2 = colorededges(P)
    graph3 = list(graph2)
    graph4 = twoBreakOnGenomeGraph(graph3, i, i1, j, j1)
    p = graphtogenome1(graph4)
    return p

def graphtogenome1(graph):
    p = list()
    cycles = findCycle1(graph)
    for nodes in cycles:
        chrom = cyclechromosome(nodes)
        p.append(chrom)
        # chrom = [str(n) if n < 0 else "+" + str(n) for n in chrom]
        # chrom = str(chrom)[1:-1].replace(",", "")
        # chrom = str(chrom)[1:-1].replace("'", "")
        # p.append("(" + chrom + ")")
    return p

def findCycle1(graph1):
    res = []
    graph = list(graph1)
    while len(graph) != 0:
        temp = []
        head = graph[0]
        graph.remove(head)
        temp.append(head)
        cur = head
        cat1 = 1
        if head[0]%2 == 0:
            cat1 = -1
        while not (cur[1] == head[0] + cat1 or cur[0] == head[0] + cat1):
            cat = 1
            if cur[1] % 2 == 0:
                cat = -1
            for g in graph:
                if g[0] == cur[1] + cat or g[1] == cur[1] + cat:
                    if g[0] == cur[1] + cat:
                        temp.append(g)
                        cur = g
                    else:
                        temp.append((g[1], g[0]))
                        cur = (g[1], g[0])
                    graph.remove(g)
                    break
        res.append(temp)
        temp = []
    res = [sum([list(r) for r in res1], []) for res1 in res]
    res = [m[-1:] + m[0:-1] for m in res]
    #print res
    return res

def ShortestRearrangementScenario(P, Q):
    print P
    rededges = list(colorededges(P))
    blueedges = list(colorededges(Q))
    while len(blueedges) != 0:
        (i, j) = blueedges[0]
        blueedges.remove((i, j))
        for r in rededges:
            if r[0] == i or r[1] == i:
                if r[0] == i:
                    i1 = r[1]
                else:
                    i1 = r[0]
                break
        for b in rededges:
            if b[0] == j or b[1] == j:
                if b[0] == j:
                    j1 = b[1]
                else:
                    j1 = b[0]
                break
        rededges = twoBreakOnGenomeGraph(rededges, i, i1, j, j1)
        p = graphtogenome1(rededges)
        #print "p", p
        p = [sorted(lst, key=abs) for lst in p]
        P = twoBreakOnGenome(p, i, i1, j, j1)
        print P


# 7.11.5
def kmer(k, s1, s2):
    sub2 = dict()
    res = []
    for i in range(0, len(s2) - k + 1):
        if not s2[i: i+k] in sub2:
            sub2[s2[i: i+k]] = []
        sub2[s2[i:i+k]].append(i)
    for j in range(0, len(s1) - k + 1):
        temp = s1[j: j+k]
        reverse = reverseStr(temp)
        if temp in sub2:
            for n in sub2[temp]:
                res.append((j, n))
        if reverse in sub2:
            for m in sub2[reverse]:
                res.append((j, m))
    return res

def reverseStr(str):
    res = ""
    for i in range(0, len(str)):
        if str[i] == 'A':
            res = res + 'T'
        elif str[i] == 'C':
            res = res + 'G'
        elif str[i] == 'G':
            res = res + 'C'
        else:
            res = res + 'A'
    res = res[::-1]
    return res

#P = [[+1, -2, -4, +3]]
#print twoBreakOnGenome(P,75, 73, 77, 80)

infile = '/Users/yijia/Documents/current work/Bio/dataset_38491_8.txt'
with open(infile) as f:
    input = f.readlines()
#k = int(input[0].strip())
#s1 = input[1].strip()
#s2 = input[2].strip()
#res = kmer(k ,s1, s2)
#num = input[0].strip()[1:-1].split(' ')
#num = [[int(n) for n in num]]
#print twoBreakOnGenome(num,75, 73, 77, 80)
#s1 = input[0].strip()[1:-1].split(')(')
#s2 = input[1].strip()[1:-1].split(')(')
#p = [[int(r) for r in s.split(' ')] for s in s1]
#q = [[int(r) for r in s.split(' ')] for s in s2]
#s1 = sum(p, [])
#s2 = sum(q, [])
#s1 = [int(n) for n in p]
#s2 = [int(n) for n in q]
#nums = [int(n) for n in nums]
#p = [[int(r) for r in s.split(', ')] for s in nums]
#pp = []
#for k in p:
#   pp = pp + k
#res = graphtogenome(pp)
#print res
#res = cyclechromosome(nums)
#res = twoBreakDistance(p, q)
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38491_8_res.txt', 'w')
#outfile.write(listToStr(res))
# for i in res:
#     outfile.write(str(i))


