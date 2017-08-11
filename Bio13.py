import sys
sys.setrecursionlimit(100000000) # 10000 is an example, try with different values

# 13.1.5
def findBWT(text):
    words = []
    for i in range(0, len(text)):
        words.append(text[i:] + text[0:i])
    words.sort()
    res = ""
    for w in words:
        res += w[-1]
    return res

#13.3.10
def transformBWT(text):
    first = sorted(text)
    last = list(text)
    preprocess(first)
    preprocess(last)
    res = first[0]
    cur = first[0]
    while len(res) < len(text):
        idx = last.index(cur)
        nex = first[idx]
        res = res + nex[0]
        cur = nex
    res = res[1:] + res[0]
    return res

def preprocess(first):
    cntA = 0
    cntC = 0
    cntG = 0
    cntT = 0
    res = []
    for i in range(0, len(first)):
        if first[i] == 'A':
            res.append(first[i] + str(cntA))
            cntA += 1
        elif first[i] == 'C':
            res.append(first[i] + str(cntC))
            cntC += 1
        elif first[i] == 'G':
            res.append(first[i] + str(cntG))
            cntG += 1
        elif first[i] == 'T':
            res.append(first[i] + str(cntT))
            cntT += 1
        else:
            res.append("$")
    return res

def findIdx(last, cur, visit):
    for i in range(len(last)):
        if last[i] == cur and visit[i] == 0:
            visit[i] = 1
            return i
    return len(visit) + 1

# 13.4.8
def BWMatching(text, pattern):
    res = []
    # print len(text)
    first = sorted(text)
    last = list(text)
    first_T = preprocess(first)
    last_T = preprocess(last)
    # print first_T
    # print first
    # print len(first_T), len(first)
    LtoF = [first_T.index(last_T[i]) for i in range(0, len(last))]
    for p in pattern:
        res.append(BWMatch(first, last, list(p), LtoF))
    return res


def BWMatch(first, last, pattern, lastToFirst):
    # print lastToFirst
    top = 0
    bottom = len(last) - 1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern.pop()
            if symbol in last[top: bottom+1]:
                topIdx = findFirst(last, top, bottom, symbol)
                bottomIdx = findLast(last, top, bottom, symbol)
                top = lastToFirst[topIdx]
                bottom = lastToFirst[bottomIdx]
            else:
                return 0
        else:
            return bottom - top + 1


def findFirst(arr, top, bottom, symbol):
    arr_copy = list(arr)
    idx = len(arr) + 1
    for i in range(top, bottom+1):
        if arr_copy[i] == symbol:
            idx = i
            return idx

def findLast(arr, top, bottom, symbol):
    arr_copy = list(arr)
    idx = len(arr) + 1
    for i in range(bottom, top-1, -1):
        if arr_copy[i] == symbol:
            idx = i
            return idx

# 13.5.7
def BetterMatching(text, pattern):
    res = []
    first = firstOccur(text)
    last = list(text)
    # first_T = preprocess(first)
    # last_T = preprocess(last)
    # LtoF = [first_T.index(last_T[i]) for i in range(0, len(last))]
    count = getCount(text)
    for p in pattern:
        res.append(BetterBWMatch(first, last, list(p), count))
    return res

def firstOccur(s):
    first = sorted(s)
    res = [first.index(s[i]) for i in range(0, len(s))]
    return first

def getCount(text):
    count = dict()
    list_t = list(text)
    for t in set(text):
        count[t] = [list_t[0:i].count(t) for i in range(0, len(text) + 1)]
    return count

def BetterBWMatch(first, last, pattern, count):
    top = 0
    bottom = len(last) - 1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern.pop()
            if symbol in last[top: bottom + 1]:
                top = first.index(symbol) + count[symbol][top]
                bottom = first.index(symbol) + count[symbol][bottom+1] - 1
            else:
                return 0
        else:
            return bottom - top + 1

# 13.7.4
def BetterMatchingIndex(text, pattern):
    checkpoint = 100
    res = []
    first = firstOccur(text)
    last = list(text)
    count = getCountPoint(text, checkpoint)
    for p in pattern:
        BetterBWMatchIdx(first, last, list(p), count, res, checkpoint)
    return sorted(res)

def getCountPoint(text, checkpoint):
    count = dict()
    list_t = list(text)
    for t in set(text):
        count[t] = [list_t[0:i].count(t) for i in range(0, len(text) + 1) if i%checkpoint == 0]
    return count

def BetterBWMatchIdx(first, last, pattern, count, res, checkpoint):
    top = 0
    bottom = len(last) - 1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern.pop()
            if symbol in last[top: bottom + 1]:
                top = first.index(symbol) + findCnt(count, symbol, top, checkpoint, last)
                bottom = first.index(symbol) + findCnt(count, symbol, top, checkpoint) - 1
        else:
            res.append(top)

def findCnt(count, symbol, top, checkpoint, last):
    if top%checkpoint == 0:
        return count[symbol][top]
    res = count[symbol][top/checkpoint]
    for i in range(top, top/checkpoint+1):
        if last[i] == symbol:
            res += 1
    return res

# 13.8.10
def matchSub(text, pattern, mis):
    res = []
    for i in range(0, len(text)):
        for p in pattern:
            if humming(str(text[i: len(p) + i]), p) <= mis:
                res.append(i)
    return sorted(res)

def humming(s1, s2):
    cnt = 0
    for i in range(0, min(len(s1), len(s2))):
        if s1[i] != s2[i]:
            cnt += 1
    cnt += max(len(s1), len(s2)) - min(len(s1), len(s2))
    return cnt

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
s = input[0].strip()
pattern = input[1].split(" ")
d = int(input[-1])
print pattern
res = matchSub(s, pattern, d)
print res
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
for r in res:
    outfile.write(str(r) + " ")
