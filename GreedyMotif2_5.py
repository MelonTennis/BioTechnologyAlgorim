import sys
import MedianString2_4
import MotifEnum2_1

# calculate probability for str and profile matrix
def pr(str, profile):
    res = 1
    for i in range(0, len(str)):
        if str[i] == 'A':
            res *= float(profile[0][i])
        elif str[i] == 'C':
            res *= float(profile[1][i])
        elif str[i] == 'G':
            res *= float(profile[2][i])
        elif str[i] == 'T':
            res *= float(profile[3][i])
    #print res
    return res

# return Profile-most Probable kmer string
def probableK(text, k, profile):
    dist = -1
    res = ""
    kmer = MedianString2_4.generateTextPattern(text, k)
    #print kmer
    for str in kmer:
        if pr(str, profile) > dist:
            dist = pr(str, profile)
            res = str
    return res

# update profile given k and motif matrix
def updateProfile(last, k):
    percent = [[0.0 for x in range(0, k)] for y in range(0, 4)]
    for line in last:
        for i in range(0, k):
            if line[i] == 'A':
                x = percent[0][i]
                percent[0][i] = 1 + x
            if line[i] == 'C':
                x = percent[1][i]
                percent[1][i] = 1 + x
            if line[i] == 'G':
                x = percent[2][i]
                percent[2][i] = 1 + x
            if line[i] == 'T':
                x = percent[3][i]
                percent[3][i] = 1 + x
    for i in range(0, len(percent)):
        for j in range(0, k):
            percent[i][j] = percent[i][j]/len(last)
    return percent

# update profile given k and motif matrix
def laplaceProfile(last, k):
    percent = [[0.0 for x in range(0, k)] for y in range(0 , 4)]
    for line in last:
        for i in range(0, k):
            if line[i] == 'A':
                x = percent[0][i]
                percent[0][i] = 1 + x
            if line[i] == 'C':
                x = percent[1][i]
                percent[1][i] = 1 + x
            if line[i] == 'G':
                x = percent[2][i]
                percent[2][i] = 1 + x
            if line[i] == 'T':
                x = percent[3][i]
                percent[3][i] = 1 + x
    for i in range(0, len(percent)):
        for j in range(0, k):
            percent[i][j] = (percent[i][j] + 1)/(len(last) + 4)
    return percent

# calculate score of motif
def score(motifs):
    res = 0
    k = len(motifs[0])
    profile = updateProfile(motifs, k)
    #print profile
    percent = [0 for x in range(0, k)]
    sequence = [0 for x in range(0, k)]
    for i in range(0, len(profile)):
        for j in range(0, k):
            if profile[i][j] > percent[j]:
                percent[j] = profile[i][j]
                sequence[j] = i
    letter = ['A', 'C', 'G', 'T']
    sequence = [letter[x] for x in sequence]
    for i in range(0, len(motifs)):
        for j in range(0, k):
            if motifs[i][j] != sequence[j]:
                res += 1
    return res

def GreedyMotif(Dna, k, t):
    best = [y[0:k] for y in Dna]
    profile = [[0.0 for x in range(0, k)] for x in range(0, 4)]
    for j in range(0, len(Dna[0]) - k + 1):
        motifs = []
        motifs.append(Dna[0][j: j+k])
        for i in range(1, t):
            #profile = updateProfile(motifs, k)
            # Laplace greedy method
            profile = laplaceProfile(motifs, k)
            temp = probableK(Dna[i], k, profile)
            motifs.append(temp)
        if score(motifs) < score(best):
            best = motifs
    return best

'''
infile = '/Users/yijia/Documents/current work/Bio/dataset_38406_9.txt'
with open(infile) as f:
    Dna = f.readlines()
k = int(Dna[0].strip().split(' ')[0])
t = int(Dna[0].strip().split(' ')[1])
Dna = [x.strip() for x in Dna[1:]]
motifs = GreedyMotif(Dna, k, t)
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38406_9_res.txt', 'w')
for p in motifs:
    outfile.write(p)
    outfile.write(' ')
outfile.close()
'''
