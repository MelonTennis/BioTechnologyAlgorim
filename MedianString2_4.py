import sys
import MotifEnum2_1


# return patterns in a text line
def generateTextPattern(text, k):
    pattern = []
    for j in range(0, len(text) - k + 1):
        t = text[j: j + k]
        pattern.append(t)
    pattern = list(set(pattern))
    return pattern

# d(pattern, text) = min HD(pattern, patterns in text)
def dist(pattern, text):
    res = sys.maxsize
    for x in generateTextPattern(text, len(pattern)):
        if MotifEnum2_1.Hamming(x, pattern) < res:
            res = MotifEnum2_1.Hamming(x, pattern)
    if(res == sys.maxsize):
        print "pattern =", pattern, "text =", text
    return res

# d(pattern, dna) = sum of (pattern, dnai in dna)
def distDna(pattern, Dna):
    res = 0
    for i in range(0, len(Dna)):
        if(len(Dna[i]) == 0):
            continue
        res += dist(pattern, Dna[i])
    return res

# Return motif finding results by median string
def MedianString(Dna, k):
    res = sys.maxsize
    Median = ""
    for str in sorted(MotifEnum2_1.neighbours(['A']*k, k)):
        if res > distDna(str, Dna):
            res = distDna(str, Dna)
            Median = str
    return Median


'''
infile = '/Users/yijia/Documents/current work/Bio/profile_most_probable.txt'
with open(infile) as f:
    Dna = f.readlines()
K = int(Dna[0].strip().split(' ')[0])
Dna = [x.strip() for x in Dna[1:]]
#print K, Dna
print MedianString(Dna, K)
#print patterns
outfile = open('/Users/yijia/Documents/current work/Bio/median_string_res.txt', 'w')
outfile.close()
'''