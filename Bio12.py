import sys
sys.setrecursionlimit(100000000) # 10000 is an example, try with different values

# 12.3.4
class Node(object):
    def __init__(self, l, str):
        self.str = str
        self.label = l
        self.child = [None for i in range(0, 28)]
        self.pos = -1
        self.childNum = 0
        self.color = -1

def TrieConstruction(patterns):
    cnt = 0
    root = Node(0, None)
    for pattern in patterns:
        cur = root
        for i in range(0, len(pattern)):
            if not cur.child[ord(pattern[i]) - ord('A')]:
                cnt += 1
                cur.child[ord(pattern[i]) - ord('A')] = Node(cnt, None)
            cur = cur.child[ord(pattern[i]) - ord('A')]
        cur.str = pattern
    return root

def printTrie(root, file, pre, pos):
    if not root:
        return
    if root.label != pre:
        file.write(str(pre) + "->" + str(root.label) + ":" + chr(ord('A')+pos))
    file.write("\n")
    for i in range(0, 26):
        if root.child[i]:
            printTrie(root.child[i], file, root.label, i)

# 12.3.8
def trieMatching(text, pattern):
    root = TrieConstruction(pattern)
    res = []
    for i in range(0, len(text)):
        t = text[i:]
        if find(t, root):
            res.append(i)
    return res

def find(p, node):
    if not node:
        return False
    cnt = node
    idx = 0
    while cnt:
        if idx >= len(p):
            return True
        if cnt.str:
            return True
        if cnt.child[ord(p[idx]) - ord('A')]:
            cnt = cnt.child[ord(p[idx]) - ord('A')]
            idx += 1
        else:
            break
    return False

# 12.5.4
def suffixTrieConsturction(text):
    patterns = [text[i:] for i in range(0, len(text)-1)]
    root = Node(0, None)
    cnt = 0
    for pattern in patterns:
        cur = root
        for i in range(0, len(pattern)):
            idx = ord(pattern[i]) - ord('A')
            if idx < 0 or idx > 25:
                idx = 26
            if not cur.child[idx]:
                cur.childNum += 1
                cur.child[idx] = Node(cnt, None)
            cur = cur.child[idx]
        cur.str = pattern
        cur.pos = cnt
        cnt += 1
    return root

def printSuffixTree(root, prefix):
    if not root:
        return
    if root.childNum <= 1:
        if root.childNum == 0:
            print prefix
            suffixRes.append(prefix)
        for i in range(0, 27):
            if root.child[i]:
                if i == 26:
                    p = '$'
                else:
                    p = chr(ord('A') + i)
                printSuffixTree(root.child[i], prefix + str(p))
    else:
        print prefix
        suffixRes.append(prefix)
        for i in range(0, 27):
            if root.child[i]:
                if i == 26:
                    p = '$'
                else:
                    p = chr(ord('A') + i)
                printSuffixTree(root.child[i], str(p))

# 12.5.5
def longestRepeat(text):
    global maxLength
    text = text + "$"
    root = suffixTrieConsturction(text)
    longestRepeatNode(root, "")
    return maxLength

def longestRepeatNode(node, pre):
    global maxLength
    if not node:
        return
    if node.childNum > 1:
        if len(pre) > len(maxLength):
            maxLength = str(pre)
        for i in range(0, 27):
            if node.child[i]:
                longestRepeatNode(node.child[i], pre + chr(ord('A') + i))
    else:
        for i in range(0, 27):
            if node.child[i]:
                longestRepeatNode(node.child[i], pre + chr(ord('A') + i))

# 12.5.6
sub = ""
def longestSubstring(text1, text2):
    global sub
    text = text1 + "#" + text2 + "$"
    pviot = len(text1) + 1
    # construct tree
    patterns = [text[i:] for i in range(0, len(text) - 1)]
    root = Node(0, None)
    root.color = 3
    cnt = 0
    for pattern in patterns:
        cur = root
        for i in range(0, len(pattern)):
            idx = ord(pattern[i]) - ord('A')
            if idx == -29:
                idx = 26
            elif idx == -30:
                idx = 27
            if not cur.child[idx]:
                cur.childNum += 1
                cur.child[idx] = Node(cnt, None)
            cur = cur.child[idx]
            if cnt < pviot:
                if cur.color < 2:
                    cur.color = 1
                else:
                    cur.color = 3
            else:
                if cur.color == 2 or cur.color == -1:
                    cur.color = 2
                else:
                    cur.color = 3
        cur.str = pattern
        cur.pos = cnt
        cnt += 1
    sub = ""
    longestShare(root, "")
    return sub

def longestShare(node, pre):
    global sub
    if not node:
        return
    if node.childNum > 1 and node.color > 2:
        if len(pre) > len(sub):
            sub = str(pre)
        for i in range(0, 28):
            if node.child[i]:
                if i == 26:
                    p = '$'
                elif i == 27:
                    p = '#'
                else:
                    p = chr(ord('A') + i)
                longestShare(node.child[i], pre + p)
    elif node.childNum == 1 and node.color > 2:
        for i in range(0, 28):
            if node.child[i]:
                if i == 26:
                    p = '$'
                elif i == 27:
                    p = '#'
                else:
                    p = chr(ord('A') + i)
                longestShare(node.child[i], pre + p)
    else:
        return

# 12.5.7
sub = ""
def shortestSubstring(text1, text2):
    global sub
    sub = text1
    text = text1 + "#" + text2 + "$"
    pviot = len(text1) + 1
    # construct tree
    patterns = [text[i:] for i in range(0, len(text) - 1)]
    root = Node(0, None)
    root.color = 3
    cnt = 0
    for pattern in patterns:
        cur = root
        for i in range(0, len(pattern)):
            idx = ord(pattern[i]) - ord('A')
            if idx == -29:
                idx = 26
            elif idx == -30:
                idx = 27
            if not cur.child[idx]:
                cur.childNum += 1
                cur.child[idx] = Node(cnt, None)
            cur = cur.child[idx]
            if cnt < pviot:
                if cur.color < 2:
                    cur.color = 1
                else:
                    cur.color = 3
            else:
                if cur.color == 2 or cur.color == -1:
                    cur.color = 2
                else:
                    cur.color = 3
        cur.str = pattern
        cur.pos = cnt
        cnt += 1
    shortestNonShare(root, "")
    return sub

def shortestNonShare(node, pre):
    global sub
    if not node:
        return
    if node.color < 3:
        if 1 < len(pre) < len(sub):
            print sub
            sub = str(pre)
        for i in range(0, 28):
            if node.child[i]:
                if i == 26:
                    p = '$'
                elif i == 27:
                    p = '#'
                else:
                    p = chr(ord('A') + i)
                shortestNonShare(node.child[i], pre + p)
    else:
        for i in range(0, 28):
            if node.child[i]:
                if i == 26:
                    p = '$'
                elif i == 27:
                    p = '#'
                else:
                    p = chr(ord('A') + i)
                shortestNonShare(node.child[i], pre + p)

# 12.6.2
def suffixArray(text):
    print len(text)
    prefix = [text[i:] for i in range(0, len(text))]
    sortPrefix = sorted(prefix)
    res = [prefix.index(sortPrefix[i]) for i in range(0, len(prefix))]
    return res

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
# patterns = [r.strip() for r in input[0:]]
# res = TrieConstruction(patterns)
# res = trieMatching(patterns[0], patterns[1:])
# root = suffixTrieConsturction(input[0].strip())
# suffixRes = []
# maxLength = ""
# print longestRepeat(input[0].strip())
# printSuffixTree(root, "")
# res = shortestSubstring(patterns[0], patterns[1])
# print res
print suffixArray(input[0].strip())
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
# printTrie(res, outfile, 0, -1)
# for r in suffixRes:
#     outfile.write(str(r) + "\n")