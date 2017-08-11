import sys
# mass table
mass = dict()
mass["G"] = 57; mass["A"] = 71; mass["S"] = 87; mass["P"] = 97; mass["T"] = 101; mass["V"] = 99;
mass["C"] = 103; mass["I"] = 113; mass["L"] = 113; mass["N"] = 114; mass["D"] = 115; mass["K"] = 128;
mass["Q"] = 128; mass["E"] = 129; mass["M"] = 131; mass["H"] = 137; mass["F"] = 147; mass["R"] = 156;
mass["Y"] = 163; mass["W"] = 186; mass["X"] = 4; mass["Z"] = 5
mass_value = mass.values()
mass_key = mass.keys()

#8.3.5
def graphSpectrum(spectrum):
    #res = []
    res = dict()
    spectrum = [0] + spectrum
    for s in spectrum:
        for m in spectrum:
            if m - s in mass_value:
                #path = str(s) + "->" + str(m) + ":" + mass_key[mass_value.index(m-s)]
                path = (s, m, mass_key[mass_value.index(m-s)])
                if s not in res.keys():
                    res[s] = path
    return res

# 8.3.7
def decodeSpectrum(spectrum):
    spectrum = [0] + spectrum
    graph = graphSpectrum(spectrum)
    head = graph[spectrum[0]][0]
    cur = graph[spectrum[0]][1]
    peptide = graph[spectrum[0]][2]
    while(cur != spectrum[-1]):
        peptide = peptide + graph[cur][2]
        head = graph[cur][0]
        cur = graph[cur][1]
    #peptide = peptide + graph[head][2]
    return peptide

# spectrum = [57, 71, 154, 185, 301, 332, 415, 429, 486]
# p = decodeSpectrum(spectrum)

# 8.5.5
def peptideVector(acio):
    prefix = []
    for i in range(0, len(acio)):
        prefix.append(getMass(acio[0:i+1]))
    #print prefix
    res = [1 if i in prefix else 0 for i in range(1, prefix[-1]+1)]
    return res

def getMass(string):
    res = 0
    for s in string:
        res = res + mass[s]
    return res

# 8.5.6
def vectorPeptide(vector):
    prefix = []
    for i in range(0, len(vector)):
        if vector[i] == 1:
            prefix.append(i+1)
    res = decodeSpectrum(prefix)
    return res

# 8.5.13
def peptideSequence(sequence):
    sequence = [0] + sequence
    #print sequence
    res = ""
    path = [-1 for i in range(0, len(sequence))]
    dp = [-sys.maxsize0 for i in range(0, len(sequence))]
    dp[0] = sequence[0]
    for i in range(1, len(sequence)):
        for j in mass_value:
            if i - j >= 0:
                dp[i] = max(dp[i], sequence[i] + dp[i -j])
                if sequence[i] + dp[i -j] == dp[i]:
                    path[i] = i-j
    cur = path[-1]
    idx = len(path) - 1
    while cur >= min(mass_value):
        res = res + mass_key[mass_value.index(idx - cur)]
        idx = cur
        cur = path[cur]
    res = res + mass_key[mass_value.index(idx - cur)]
    res = res[::-1]
    return res

# 8.6.2
def peptideIdentify(spectrum, proteome):
    maxV = -sys.maxsize0
    res = ""
    for ii in range(0, len(proteome)):
        for jj in range(ii, len(proteome)):
            prefix = getMass(proteome[ii: jj])
            if prefix > len(spectrum):
                break
            if prefix == len(spectrum):
                vector = peptideVector(proteome[ii: jj])
                score = sum([vector[i] * spectrum[i] for i in range(0, len(spectrum))])
                if score > maxV:
                    maxV = score
                    res = proteome[ii:jj]
    return res

# 8.6.7
def PSMsearch(spectral, proteome, threshold):
    PSM = set()
    for vector in spectral:
        peptide = peptideIdentify(vector, proteome)
        if getScore(peptide, vector) >= threshold:
            PSM.add(peptide)
    return PSM

def getScore(peptide, spectrum):
    # print peptide
    # vector = peptideVector(peptide)
    # if len(vector) != len(spectrum):
    #     return -1000000
    # else:
    #     return sum([vector[i] * spectrum[i] for i in range(0, len(spectrum))])
    prefix = getMass(peptide)
    if prefix != len(spectrum):
        return -1
    else:
        vector = peptideVector(peptide)
        score = sum([vector[i] * spectrum[i] for i in range(0, len(spectrum))])
        return score

# 8.8.3
def spectralDection(spectral, threshold, max_score):
    spectral = [0] + spectral
    table = [[0 for i in range(0, len(spectral))] for j in range(0, max_score)]
    table[0][0] = 1
    for i in range(0, len(spectral)):
        for t in range(0, max_score):
            for a in mass_value:
                if i - a >= 0 and t - spectral[i] >= 0 and t - spectral[i] < max_score:
                    #print i, t, i - a, t - spectral[i]
                    table[t][i] += table[t - spectral[i]][i - a]
    res = 0
    for i in range(threshold, max_score):
        res += table[i][len(spectral) - 1]
    return res

# 8.8.8
def probability(spectral, threshold, max_score):
    spectral = [0] + spectral
    table = [[0 for i in range(0, len(spectral))] for j in range(0, max_score)]
    table[0][0] = 1.0
    for i in range(0, len(spectral)):
        for t in range(0, max_score):
            for a in mass_value:
                if i - a >= 0 and t - spectral[i] >= 0 and t - spectral[i] < max_score:
                    #print i, t, i - a, t - spectral[i]
                    table[t][i] += float(table[t - spectral[i]][i - a])/20.0
    res = 0
    for i in range(threshold, max_score):
        res += table[i][len(spectral) - 1]
    #print table
    return res

# 8.11.3
def spectralAlignment(peptide, spectral, k):
    dp = [[[0 for i in range(0, len(spectral))] for j in range(0, len(peptide) + 1)] for kk in range(0, k+1)]
    prefix = [getMass(peptide[0:i]) for i in range(0, len(peptide))]
    print prefix
    dp[0][0][0] = 0
    for m in range(1, k+1):
        dp[m][0][0] = -sys.maxsize0
    for t in range(0, k+1):
        for i in range(0, len(peptide)+1):
            for j in range(0, len(spectral)):
                if i == 0 and j == 0 and t == 0:
                    continue
                if i != 0:
                    diff = mass[peptide[i-1]]
                else:
                    diff = 0
                #print diff
                for j1 in range(0, j):
                    if j - diff >= 0 and prefix[i-1] - diff in prefix and t >= 1 and j - diff <= len(peptide):
                        #print t, i, j, j - diff, prefix.index(prefix[i-1] - diff), j1
                        dp[t][i][j] = spectral[j] + max(dp[t][j - diff][prefix.index(prefix[i-1] - diff)], dp[t - 1][prefix.index(prefix[i-1] - diff)][j1])
    maxValue = -100
    position = (0, 0, 0)
    for d in dp:
        for mm in d:
            print mm
        print "\n"
    for q in range(0, k+1):
        maxValue = max(maxValue, dp[q][len(peptide)][len(spectral) - 1])
        if maxValue == dp[q][len(peptide)][len(spectral) - 1]:
            position = (q, len(peptide), len(spectral) - 1)
    print position
    change = [0 for i in range(0, len(peptide))]
    i = len(peptide)
    j = len(spectral) - 1
    t = position[0]
    iter = 0
    while not (i == 0 and j == 0):
        diff = mass[peptide[i - 1]]
        for j1 in range(0, j):
            if prefix[i-1] - diff in prefix:
                if j - diff >= 0 and t >= 1 and j - diff <= len(peptide) and dp[t - 1][prefix.index(prefix[i-1] - diff)][j1] + spectral[j] == dp[t][i][j]:
                    change[-1 - iter] = diff - (prefix[j] - prefix[j-1])
                    i = prefix.index(prefix[i-1] - diff)
                    j = j1
                else:
                    if i == 0:
                        break
                    change[-1 - iter] = 0
                    j = j - diff
                    i = prefix.index(prefix[i-1] - diff)
                iter += 1
    print change

# 8.11
def alignmentSpectral(peptide, score, k):
    dp = [[[0 for i in range(0, k + 1)] for j in range(0, len(score))] for m in range(0, len(peptide) + 1)]
    bt = [[[-sys.maxsize for i in range(0, k + 1)] for j in range(0, len(score))] for m in range(0, len(peptide) + 1)]
    bt[0][0][0] = 0
    for i in range(1, len(peptide) + 1):
        for j in range(0, len(score)):
            for t in range(0, k + 1):
                for l in range(0, j):
                    if j - l == mass[peptide[i - 1]]:
                        if bt[i-1][l][t] != -sys.maxsize and (bt[i][j][t] == -sys.maxsize or dp[i][j][t] < dp[i-1][l][t] + score[j]):
                            dp[i][j][t] = dp[i-1][l][t] + score[j]
                            bt[i][j][t] = j - l
                    else:
                        if t > 0 and bt[i-1][l][t-1] != -sys.maxsize and (bt[i][j][t] == -sys.maxsize or dp[i][j][t] < dp[i-1][l][t-1] + score[j]):
                            dp[i][j][t] = dp[i-1][l][t-1] + score[j]
                            bt[i][j][t] = j-l

    max_t = -1
    max_score = 0
    for i in range(0, k+1):
        if bt[len(peptide)][len(score) - 1][i] != -sys.maxsize and max_t == -1 or max_score < dp[len(peptide)][len(score) - 1]:
            max_t = i
            max_score = dp[len(peptide)][len(score) - 1][i]

    Mass = len(score) - 1
    diffs = []
    for i in range(len(peptide), 0, -1):
        diff = bt[i][Mass][max_t]
        Mass -= diff
        diffs.append(diff)
        if diff != mass[peptide[i-1]]:
            max_t -= 1
    diffs.reverse()
    #print diffs

    res = ""
    for i in range(0, len(peptide)):
        res += peptide[i]
        dif = diffs[i] - mass[peptide[i]]
        if dif != 0:
            res += "("
            if dif > 0:
                res += "+"
            else:
                res += ""
            res += str(dif)
            res += ')'
    return res

p = "XXZ"
spect = [int(m) for m in "4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 -1".split(' ')]
k = 2
print alignmentSpectral(p, spect, k)

# spectrum = [4, -3, -2, 3, 3, -4, 5, -3, -1, -1, 3, 4, 1, 3]
# thre = 1
# maxV = 8
# print probability(spectrum, thre, maxV)

# spectral = [[-1, 5, -4, 5, 3, -1, -4, 5, -1, 0, 0, 4, -1, 0, 1, 4, 4],[-4, 2, -2, -4, 4, -5, -1, 4, -1, 2, 5, -3, -1, 3, 2, -3]]
# pro = "XXXZXZXXZXZXXXZXXZX"
# thre = 5
# print PSMsearch(spectral, pro, thre)

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
p = input[0].strip()
s1 = input[1].strip().split(" ")
s1 = [int(m) for m in s1]
k = int(input[2].strip())
#print s1, t, k
res = alignmentSpectral(p, s1, k)
print res
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
outfile.write(res)
