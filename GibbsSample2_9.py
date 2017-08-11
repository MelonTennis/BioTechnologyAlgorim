import sys
import numpy
import random
import MedianString2_4
import MotifEnum2_1
import GreedyMotif2_5
import RandomMotif2_7


def weightChoice(choices, weights):
    s = sum(weights)
    ts = random.uniform(0, s)
    for m, w in enumerate(weights):
        s -= w
        if s < ts:
            return choices[m]

def randomChoose(text, k, profile):
    prob = []
    pp = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i: i + k]
        p = GreedyMotif2_5.pr(pattern, profile)
        prob.append(p)
        pp.append(pattern)
    prob = [x/sum(prob) for x in prob]
    return weightChoice(pp, prob)

def Gibbs(Dna, k, t, N):
    motifs = RandomMotif2_7.randomChoose(Dna, k)
    best = motifs
    for j in range(0, N):
        i = random.randint(0, t-1)
        profile = GreedyMotif2_5.laplaceProfile(motifs[0: i] + motifs[i+1:], k)
        motifs[i] = randomChoose(Dna[i], k, profile)
        #print GreedyMotif2_5.score(motifs)
        if GreedyMotif2_5.score(motifs) < GreedyMotif2_5.score(best):
            #print GreedyMotif2_5.score(motifs)
            best = motifs
    return best

'''
infile = '/Users/yijia/Documents/current work/Bio/dataset_38409_4.txt'
with open(infile) as f:
    Dna = f.readlines()
k = int(Dna[0].strip().split(' ')[0])
t = int(Dna[0].strip().split(' ')[1])
N = int(Dna[0].strip().split(' ')[2])
Dna = [x.strip() for x in Dna[1:]]
motifs = Gibbs(Dna, k, t, N)
print GreedyMotif2_5.score(motifs)
for i in range(0, 20):
    tempMotifs = Gibbs(Dna, k, t, N)
    if GreedyMotif2_5.score(motifs) > GreedyMotif2_5.score(tempMotifs):
        motifs = tempMotifs
print GreedyMotif2_5.score(motifs)
infile2 = '/Users/yijia/Documents/current work/Bio/sample.txt'
with open(infile2) as f2:
    sample = f2.readlines()
sample = [x.strip() for x in sample]
print GreedyMotif2_5.score(sample)
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38409_4_res.txt', 'w')
for p in motifs:
    outfile.write(p)
    outfile.write('\n')
outfile.close()
'''