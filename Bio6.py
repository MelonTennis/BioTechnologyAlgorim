import re

def readMatrix(infile):
    with open(infile) as f:
        line = f.readlines()
    alphabet = dict()
    word = re.split('\s+', line[0].strip())
    for i in range(0, len(word)):
        alphabet[word[i]] = i
    matrix = [re.split('\s+', l.strip()) for l in line[1:]]
    dictionary = dict()
    for m in matrix:
        dictionary[m[0]] = [int(x) for x in m[1:]]
    return (alphabet, dictionary)

#infile = '/Users/yijia/Documents/current work/Bio/BLOSUM62.txt'
#(a, m) = readMatrix(infile)

# get point from BLOSUM62 score matrix
def getPoint(a, m, c1, c2):
    index = a[c2]
    #print c1, c2, m[c1][index]
    return m[c1][index]

# 6.2.3
def globalAlignment(s1, s2, alphabet, matrix):
    penalty = 5
    #if len(s2) > len(s1):
    #    return globalAlignment(s2, s1, alphabet, matrix)
    dp = [[0 for i in range(0, len(s2)+1)] for j in range(0, len(s1)+1)]
    for i in range(0, len(s1)+1):
        if i == 0:
            for j in range(1, len(s2)+1):
                dp[i][j] = dp[i][j-1] - penalty
        else:
            for j in range(0, len(s2)+1):
                if j == 0:
                    dp[i][j] = dp[i-1][j] - penalty
                else:
                    dp[i][j] = max(dp[i-1][j] - penalty, dp[i][j-1] - penalty)
                    point = getPoint(alphabet, matrix, s1[i - 1], s2[j - 1])
                    dp[i][j] = max(dp[i][j], dp[i-1][j-1] + point)
    (align1, align2) = constructAlign(alphabet, matrix, dp, s1, s2, penalty)
    #print dp
    #print dp[len(s1)][len(s2)]
    return (dp[len(s1)][len(s2)], align1, align2)

# 6.2.3 construct alignment
def constructAlign(alphabet, matrix, dp, s1, s2, penalty):
    align1 = ""
    align2 = ""
    i = len(dp)-1
    j = len(dp[0])-1
    while i != 0 and j != 0 and dp[i][j] != 0:
        if dp[i][j] == dp[i-1][j] - penalty:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i = i-1
        elif dp[i][j] == dp[i][j-1] - penalty:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j = j-1
        else:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i = i-1
            j = j-1
    if i != 0 and j == 0:
        align1 = s1[0] + align1
        align2 = '-' + align2
    elif i == 0 and j != 0:
        align1 = '-' + align1
        align2 = s2[0] + align2
    return (align1, align2)

#infileL = '/Users/yijia/Documents/current work/Bio/PAM250.txt'
#(al, ml) = readMatrix(infileL)

# 6.2.10 local alignment problem
def localAlign(s1, s2, alphabet, matrix):
    penalty = 5
    mx = 0
    my = 0
    maxValue = -100000
    dp = [[0 for i in range(0, len(s2) + 1)] for j in range(0, len(s1) + 1)]
    for i in range(0, len(s1) + 1):
        if i == 0:
            for j in range(1, len(s2) + 1):
                dp[i][j] = max(dp[i][j - 1] - penalty, 0)
        else:
            for j in range(0, len(s2) + 1):
                if j == 0:
                    dp[i][j] = max(dp[i - 1][j] - penalty, 0)
                else:
                    dp[i][j] = max(dp[i - 1][j] - penalty, dp[i][j - 1] - penalty, 0)
                    point = getPoint(alphabet, matrix, s1[i - 1], s2[j - 1])
                    dp[i][j] = max(dp[i][j], dp[i - 1][j - 1] + point, 0)
        maxValue = max(maxValue, max(dp[i]))
        if(maxValue == max(dp[i])):
            my = dp[i].index(maxValue)
            mx = i
    (align1, align2) = constructAlign(alphabet, matrix, dp, s1, s2, penalty, (mx, my))
    #print dp
    #print maxValue
    return (maxValue, align1, align2)

def constructAlign(alphabet, matrix, dp, s1, s2, penalty, (mx, my)):
    align1 = ""
    align2 = ""
    i = mx
    j = my
    while i != 0 and j != 0 and dp[i][j] != 0:
        if dp[i][j] == dp[i-1][j] - penalty:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i = i-1
        elif dp[i][j] == dp[i][j-1] - penalty:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j = j-1
        else:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i = i-1
            j = j-1
    if i != 0 and j == 0:
        align1 = s1[0] + align1
        align2 = '-' + align2
    elif i == 0 and j != 0:
        align1 = '-' + align1
        align2 = s2[0] + align2
    return (align1, align2)

# 6.3.3 Edit Distance
def editDistance(s1, s2):
    dp = [[0 for i in range(0, len(s2)+1)] for j in range(0, len(s1)+1)]
    for ii in range(0, len(s1)+1):
        dp[ii][0] = ii
    for jj in range(0, len(s2)+1):
        dp[0][jj] = jj
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i][j - 1] + 1, dp[i - 1][j] + 1, dp[i - 1][j - 1] + cost)
    #(score, align1, align2) = constructAlign(dp, s1, s2)
    return dp[len(s1)][len(s2)]

def fittingAlignment(s1, s2):
    mx = 0
    maxValue = -110000
    dp = [[0 for i in range(0, len(s2) + 1)] for j in range(0, len(s1) + 1)]
    for i in range(0, len(s1) + 1):
        if i == 0:
            for j in range(1, len(s2) + 1):
                dp[i][j] = dp[i][j - 1] -1
        else:
            for j in range(0, len(s2) + 1):
                if j == 0:
                    dp[i][j] = max(dp[i - 1][j] - 1, 0)
                else:
                    dp[i][j] = max(dp[i - 1][j] - 1, dp[i][j - 1] - 1)
                    point = 1 if s1[i-1] == s2[j-1] else -1
                    dp[i][j] = max(dp[i][j], dp[i - 1][j - 1] + point)
                if j == len(s2) and i == len(s1):
                    dp[len(s1)][len(s2)] = max([x[len(s2)] for x in dp])
    for k in range(0, len(s1)+1):
        if dp[k][len(s2)] > maxValue:
            maxValue = dp[k][len(s2)]
            mx = k
    print maxValue
    (align1, align2) = constructAlign(dp, s1, s2, mx, len(s2))
    return (maxValue, align1, align2)

def constructAlign(dp, s1, s2, mx, my, penalty):
    align1 = ""
    align2 = ""
    i = mx
    j = my
    while j != 0:
        if dp[i][j] == dp[i-1][j-1] - penalty:
            align1 = s1[i - 1] + align1
            align2 = s2[j - 1] + align2
            i = i - 1
            j = j - 1
        elif dp[i][j] == dp[i-1][j] - penalty:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i = i-1
        elif dp[i][j] == dp[i][j-1] - penalty:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j = j-1
        else:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i = i-1
            j = j-1
    return (align1, align2)

def overlapAlign(s1, s2):
    my = 0
    maxValue = -110000
    penalty = 2
    dp = [[0 for i in range(0, len(s2) + 1)] for j in range(0, len(s1) + 1)]
    for i in range(0, len(s1) + 1):
        if i == 0:
            for j in range(1, len(s2) + 1):
                dp[i][j] = dp[i][j - 1] - penalty
        else:
            for j in range(0, len(s2) + 1):
                if j == 0:
                    dp[i][j] = max(dp[i - 1][j] - penalty, 0)
                else:
                    dp[i][j] = max(dp[i - 1][j] - penalty, dp[i][j - 1] - penalty)
                    point = 1 if s1[i - 1] == s2[j - 1] else -penalty
                    dp[i][j] = max(dp[i][j], dp[i - 1][j - 1] + point)
                if i == len(s1) and j == len(s2):
                    dp[len(s1)][len(s2)] = max(dp[len(s1)])
    for k in range(0, len(s2) + 1):
        if dp[len(s1)][k] > maxValue:
            maxValue = dp[len(s1)][k]
            my = k
    #for d in dp:
    #    print d
    print maxValue
    (align1, align2) = constructAlign(dp, s1, s2, len(s1), my, penalty)
    #print align1, align2
    return (maxValue, align1, align2)

# 6.4.8
def affineGap(s1, s2, alphabet, matrix):
    lower = [[0 for i in range(0, len(s2)+1)] for j in range(0, len(s1)+1)]
    upper = [[0 for i in range(0, len(s2)+1)] for j in range(0, len(s1)+1)]
    middle = [[0 for i in range(0, len(s2)+1)] for j in range(0, len(s1)+1)]
    for i in range(0, len(s1)+1):
        upper[i][0] = -1000000
    for i in range(0, len(s2)+1):
        lower[0][i] = -1000000
    gapPenalty = 11
    extensionPenalty = 1
    for i in range(0, len(s1)+1):
        for j in range(0, len(s2)+1):
            if i == 0 and j == 0:
                continue
            if i == 0:
                upper[i][j] = max(upper[i][j-1] - extensionPenalty, middle[i][j-1] - gapPenalty)
                middle[i][j] = max(lower[i][j], upper[i][j])
            else:
                if j == 0:
                    lower[i][j] = max(lower[i-1][j] - extensionPenalty, middle[i-1][j] - gapPenalty)
                    middle[i][j] = max(lower[i][j], upper[i][j])
                else:
                    lower[i][j] = max(lower[i-1][j] - extensionPenalty, middle[i-1][j] - gapPenalty)
                    upper[i][j] = max(upper[i][j-1] - extensionPenalty, middle[i][j-1] - gapPenalty)
                    point = getPoint(alphabet, matrix, s1[i - 1], s2[j - 1])
                    middle[i][j] = max(lower[i][j], middle[i-1][j-1] + point, upper[i][j])
    print middle[len(s1)][len(s2)]
    (align1, align2) = constructAlign(lower, middle, upper, s1, s2, extensionPenalty, gapPenalty)
    #print align1, align2
    return (middle[len(s1)][len(s2)], align1, align2)

def constructAlign(lower, middle, upper, s1, s2, extension, gap):
    align1 = ""
    align2 = ""
    i = len(s1)
    j = len(s2)
    while i != 0 and j != 0:
        if middle[i][j] == lower[i][j]:
            (align1, align2, i, j) = findLower(lower, middle, align1, align2, s1, s2, extension, gap, i, j)
        elif middle[i][j] == upper[i][j]:
            (align1, align2, i, j) = findUpper(upper, middle, align1, align2, s1, s2, extension, gap, i, j)
        else:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i = i - 1
            j = j - 1
    if i != 0 and j == 0:
        align1 = s1[0] + align1
        align2 = '-' + align2
    elif i == 0 and j != 0:
        align1 = '-' + align1
        align2 = s2[0] + align2
    return (align1, align2)

def findLower(lower, middle, align1, align2, s1, s2, extension, gap, i, j):
    while lower[i][j] == lower[i-1][j] - extension:
        align1 = s1[i-1] + align1
        align2 = '-' + align2
        i = i-1
    align1 = s1[i-1] + align1
    i = i - 1
    align2 = '-' + align2
    return (align1, align2, i, j)

def findUpper(upper, middle, align1, align2, s1, s2, extension, gap, i, j):
    while upper[i][j] == upper[i][j-1] - extension:
        align2 = s2[j-1] + align2
        align1 = '-' + align1
        j = j-1
    align2 = s2[j-1] + align2
    j = j - 1
    align1 = '-' + align1
    return (align1, align2, i, j)

def middleEdge(s1, s2, alphabet, matrix):
    middleCol = len(s2)/2
    penalty = 5
    middleList1 = fromSource(s1, s2, middleCol,alphabet, matrix, penalty)
    (middleList2, nextplace) = fromSink(s1, s2, middleCol, alphabet, matrix, penalty)
    middleList = [middleList1[i] + middleList2[i] for i in range(0, len(middleList2))]
    index = middleList.index(max(middleList))
    middlePoint = (index, middleCol)
    nextPoint = (len(s1) - nextplace[index][0], len(s2) - nextplace[index][1])
    return (middlePoint, nextPoint)

def fromSource(s1, s2, middleCol, alphabet, matrix, penalty):
    dp = [[0 for i in range(0, middleCol+1)] for j in range(0, len(s1)+1)]
    for i in range(0, len(s1) + 1):
        if i == 0:
            for j in range(1, middleCol + 1):
                dp[i][j] = dp[i][j - 1] - penalty
        else:
            for j in range(0, middleCol + 1):
                if j == 0:
                    dp[i][j] = dp[i - 1][j] - penalty
                else:
                    dp[i][j] = max(dp[i - 1][j] - penalty, dp[i][j - 1] - penalty)
                    point = getPoint(alphabet, matrix, s1[i - 1], s2[j - 1])
                    dp[i][j] = max(dp[i][j], dp[i - 1][j - 1] + point)
    return [d[middleCol] for d in dp]

def fromSink(s1, s2, middleCol, alphabet, matrix, penalty):
    s1_reverse = s1[::-1]
    s2_reverse = s2[::-1]
    middle = middleCol
    place = dict()
    middleCol = len(s2) - middleCol
    dp = [[0 for i in range(0, middleCol+1)] for j in range(0, len(s1)+1)]
    for i in range(0, len(s1) + 1):
        if i == 0:
            for j in range(1, middleCol + 1):
                dp[i][j] = dp[i][j - 1] - penalty
                if j == middleCol:
                    place[len(s1)] = (0, middleCol)
        else:
            for j in range(0, middleCol + 1):
                if j == 0:
                    dp[i][j] = dp[i - 1][j] - penalty
                else:
                    dp[i][j] = max(dp[i - 1][j] - penalty, dp[i][j - 1] - penalty)
                    point = getPoint(alphabet, matrix, s1[i - 1], s2[j - 1])
                    dp[i][j] = max(dp[i][j], dp[i - 1][j - 1] + point)
                    if j == middleCol:
                        if dp[i][j] == dp[i - 1][j] - penalty:
                            place[len(s1) - i] = (i-1, j)
                        elif dp[i][j] == dp[i][j - 1] - penalty:
                            place[len(s1) - i] = (i, j-1)
                        else:
                            place[len(s1) - i] = (i-1, j-1)
    #print place
    res = [d[middleCol] for d in dp]
    res.reverse()
    return (res, place)

# 6.5.13
def linearSpaceAlignment(v, w, top, bottom, left, right, alphabet, matrix):
    if left == right:
        return (bottom - top, v[top + 1: bottom + 1], '-'*(bottom - top))
    if top == bottom:
        return (right - left, '-'*(right - left), w[left + 1: right + 1])
    middle = (left + right)/2
    midNode = findMidNode(v, w, top, bottom, left, right, alphabet, matrix)
    midEdge = findMidEdge(v, w, top, bottom, left, right, alphabet, matrix)
    print linearSpaceAlignment(v, w, top, midNode, left, middle, alphabet, matrix)
    s1 = v[top+1 :bottom+1]
    s2 = w[left+1 :right+1]
    if len(s1) != 0 and len(s2) != 0:
        (mid, edge) = middleEdge(s1, s2, alphabet, matrix)
        print (s1[mid[0]-1] + s1[edge[0] -1], s2[mid[1] - 1] + s2[edge[0] - 1])
    if midEdge == 'D' or midEdge == 'H':
        middle = middle + 1
    if midEdge == 'V' or midEdge == 'D':
        midNode = midNode + 1
    print linearSpaceAlignment(v, w, midNode, bottom, middle, right, alphabet, matrix)
    return

def findMidNode(v, w, top, bottom, left, right, alphabet, matrix):
    s1 = v[top+1:bottom+1]
    s2 = w[left+1:right+1]
    if len(s1) == 0 or len(s2) == 0:
        return top
    mid = middleEdge(s1, s2, alphabet, matrix)[0][0]
    return mid + top

def findMidEdge(v, w, top, bottom, left, right, alphabet, matrix):
    s1 = v[top+1:bottom+1]
    s2 = w[left+1:right+1]
    if len(s1) == 0 or len(s2) == 0:
        return 'D'
    (mid, edge) = middleEdge(s1, s2, alphabet, matrix)
    if edge[0] == mid[0]:
        return 'H'
    elif edge[1] == mid[1]:
        return 'V'
    else:
        return 'D'

# 6.6
def multipleAlignmet(s1, s2, s3):
    s = [[[0 for i in range(0, len(s3)+1)] for j in range(0, len(s2)+1)] for w in range(0, len(s1)+1)]
    for j in range(1, len(s2) + 1):
        for k in range(1, len(s3) + 1):
            s[0][j][k] = max(s[0][j-1][k-1], s[0][j][k-1], s[0][j-1][k-1])
    for i in range(1, len(s1) + 1):
        for k in range(1, len(s3) + 1):
            s[i][0][k] = max(s[i-1][0][k], s[i][0][k-1], s[i-1][0][k-1])
    for i in range(0, len(s1) + 1):
        for j in range(0, len(s2) + 1):
            s[i][j][0] = max(s[i-1][j][0], s[i][j-1][0], s[i-1][j-1][0])
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            for k in range(1, len(s3) + 1):
                    point = 0
                    if s1[i-1] == s2[j-1] and s2[j-1] == s3[k-1]:
                        point = 1
                    #print point
                    s[i][j][k] = max(s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1], point + s[i-1][j-1][k-1])
    #print s[len(s1)][len(s2)][len(s3)]
    (align1, align2, align3) = constructAlign3D(s, s1, s2, s3)
    return (s[len(s1)][len(s2)][len(s3)], align1, align2, align3)

def constructAlign3D(s, s1, s2, s3):
    align1 = ""
    align2 = ""
    align3 = ""
    i = len(s1)
    j = len(s2)
    k = len(s3)
    while i != 0 and j != 0 and k != 0:
        if s[i][j][k] == s[i-1][j][k]:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            align3 = '-' + align3
            i = i -1
        elif s[i][j][k] == s[i][j-1][k]:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            align3 = '-' + align3
            j = j-1
        elif s[i][j][k] == s[i][j][k-1]:
            align1 = '-' + align1
            align2 = '-' + align2
            align3 = s3[k-1] + align3
            k = k-1
        elif s[i][j][k] == s[i-1][j-1][k]:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            align3 = '-' + align3
            i = i-1
            j = j-1
        elif s[i][j][k] == s[i-1][j][k-1]:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            align3 = s3[k-1] + align3
            i = i-1
            k = k-1
        elif s[i][j][k] == s[i][j-1][k-1]:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            align3 = s3[k-1] + align3
            j = j-1
            k = k-1
        else:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            align3 = s3[k-1] + align3
            i = i-1
            j = j-1
            k = k-1
    if i != 0 or j != 0 or k != 0:
        if i != 0:
            align1 = s1[i-1] + align1
        else:
            align1 = '-' + align1
        if j != 0:
            align2 = s2[j-1] + align2
        else:
            align2 = '-' + align2
        if k != 0:
            align3 = s3[k-1] + align3
        else:
            align3 = '-' + align3
   # print len(align1), len(align2), len(align3)
    return (align1, align2, align3)


#print multipleAlignmet("ATATCCG", "TCCGA", "ATGTACTG")

infile = '/Users/yijia/Documents/current work/Bio/dataset_38473_5.txt'
with open(infile) as f:
    input = f.readlines()
s1 = input[0].strip()
s2 = input[1].strip()
s3 = input[2].strip()
(length, string1, string2, string3) = multipleAlignmet(s1, s2, s3)
#print middleEdge(s1, s2, a, m)
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38473_5_res.txt', 'w')
outfile.write(str(length))
outfile.write('\n')
outfile.write(string1)
outfile.write('\n')
outfile.write(string2)
outfile.write('\n')
outfile.write(string3)
