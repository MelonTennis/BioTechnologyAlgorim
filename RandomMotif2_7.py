import sys
import random
import MedianString2_4
import MotifEnum2_1
import GreedyMotif2_5

# return a randomly chosen from each line from Dna motifs
def randomChoose(Dna, k):
    motifs = []
    for i in range(0, len(Dna)):
        rnd = random.randint(0, len(Dna[i]) - k)
        motifs.append(Dna[i][rnd: rnd+k])
    return motifs

# get profile most motifs from profile and Dna
def makeMotifs(profile, Dna):
    motifs = []
    for line in Dna:
        temp = GreedyMotif2_5.probableK(line, k, profile)
        motifs.append(temp)
    return motifs

# random motif choosing
def RandomMotif(Dna, k, t):
    motifs = randomChoose(Dna, k)
    best = motifs
    while True:
        profile = GreedyMotif2_5.laplaceProfile(motifs, k)
        motifs = makeMotifs(profile, Dna)
        if GreedyMotif2_5.score(motifs) < GreedyMotif2_5.score(best):
            best = motifs
        else:
            #print "bestscore", GreedyMotif2_5.score(best)
            return best

'''
infile = '/Users/yijia/Documents/current work/Bio/dataset_38407_5.txt'
with open(infile) as f:
    Dna = f.readlines()
k = int(Dna[0].strip().split(' ')[0])
t = int(Dna[0].strip().split(' ')[1])
Dna = [x.strip() for x in Dna[1:]]
motifs = RandomMotif(Dna, k, t)
for i in range(0, 999):
    tempMotifs = RandomMotif(Dna, k, t)
    if GreedyMotif2_5.score(motifs) > GreedyMotif2_5.score(tempMotifs):
        motifs = tempMotifs
print GreedyMotif2_5.score(motifs)
#print motifs

infile2 = '/Users/yijia/Documents/current work/Bio/sample.txt'
with open(infile2) as f2:
    sample = f2.readlines()
sample = [x.strip() for x in sample]
print GreedyMotif2_5.score(sample)

outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38407_5_res.txt', 'w')
for p in motifs:
    outfile.write(p)
    outfile.write('\n')
outfile.close()
'''