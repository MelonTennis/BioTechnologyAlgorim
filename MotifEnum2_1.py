
# generate every k-mer pattern in DNA
def GeneratePattern(Dna, k):
    pattern = []
    for i in range(0, len(Dna)):
        for j in range(0, len(Dna[i]) - k + 1):
            t = Dna[i][j: j+k]
            pattern.append(t)
    pattern = list(set(pattern))
    #print "pattern", pattern
    return pattern


# Hamming Distance from two String
def Hamming(p1, p2):
    assert(len(p1) == len(p2))
    count = 0
    for i in range(0, len(p1)):
        if p1[i] != p2[i]:
            count += 1
    #print p1, p2, count
    return count

# Return if p appears in all Dna
def appear(Dna, p, d, k):
    for i in range(0, len(Dna)):
        line = False
        temp = False
        for j in range(0, len(Dna[i]) - k + 1):
            if Hamming(Dna[i][j: j + k], p) <= d:
                temp = True
            line = line or temp
        if line == False:
            return False
    #print "p", p
    return True

# Return all d mistakes neighbours for pattern
def neighbours(pattern, d):
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']
    res = []
    suffix = neighbours(pattern[1:], d)
    for str in suffix:
        if Hamming(str, pattern[1:]) < d:
            for x in ['A', 'C', 'G', 'T']:
                res.append(x + str)
        else:
            res.append(pattern[0] + str)
    return list(set(res))


# Generate all (k, d) motifs in Dna
def MotifEnum(Dna, k, d):
    patterns = []
    kmer = GeneratePattern(Dna, k)
    for p1 in kmer:
        for x in neighbours(p1, d):
            if appear(Dna, x, d, k):
                patterns.append(x)
    patterns = list(set(patterns))
    return patterns


'''
infile = '/Users/yijia/Documents/current work/Bio/dataset_38419_3.txt'
with open(infile) as f:
    Dna = f.readlines()
K = int(Dna[0].strip().split(' ')[0])
Dna = [x.strip() for x in Dna[1:]]
print K, Dna
patterns = GeneratePattern(Dna, K)
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38419_3_res.txt', 'w')
for p in patterns:
    outfile.write(p)
    outfile.write('\n')
outfile.close()
'''