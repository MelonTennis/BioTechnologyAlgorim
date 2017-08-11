#Chapter 4 for Bioinformation
#4.2.4 Protein translation problem
def proteinTrans(Rna):
    infile = '/Users/yijia/Documents/current work/Bio/RNA_codon.txt'
    codon = dict()
    with open(infile) as f:
        RNA = f.readlines()
    Pair = [x.strip() for x in RNA]
    for pair in Pair:
        if len(pair) == 3:
            codon[pair] = ""
        else:
            codon[pair.split(' ')[0]] = pair.split(' ')[1]
    protin = ""
    i = 0
    while i < len(Rna):
        if Rna[i: i + 3] != "":
            protin = protin + codon[Rna[i: i + 3]]
        else:
            break
        i = i + 3
    return protin

# reverse complement of Dna sequence
def reverseDna(Dna):
    res = list()
    for i in range(0, len(Dna)):
        if Dna[i] == 'A':
            res.append('T')
        elif Dna[i] == 'C':
            res.append('G')
        elif Dna[i] == 'G':
            res.append("C")
        else:
            res.append("A")
    res.reverse()
    return "".join(str(x) for x in res)

# transform Dna to Rna
def DtoR(Dna):
    res = list()
    for i in range(0, len(Dna)):
        if Dna[i] == 'A':
            res.append('U')
        elif Dna[i] == 'C':
            res.append('G')
        elif Dna[i] == 'G':
            res.append("C")
        else:
            res.append("A")
    res.reverse()
    return "".join(str(x) for x in res)

# find peptide
def peptideEncode(Dna, protin):
    res = list();
    for i in range(0, len(Dna) - 3*len(protin)):
        sequence = Dna[i: i + 3*len(protin)]
        sequen_rev = reverseDna(sequence)
        ori_Rna = DtoR(sequence)
        rev_Rna = DtoR(ori_Rna)
        if proteinTrans(ori_Rna) == protin or proteinTrans(rev_Rna) == protin:
            res.append(sequence)
    return res

# calculate number of peptide 4.6.3
def subpeNum(n):
    return n*(n+1)/2+1

def spectometer(Dna):
    infile = '/Users/yijia/Documents/current work/Bio/integer_mass_table.txt'
    mass = dict()
    with open(infile) as f:
        alph = f.readlines()
    Pair = [x.strip() for x in alph]
    for pair in Pair:
        mass[pair.split(' ')[0]] = int(pair.split(' ')[1])
    res = list()
    for i in range(0, len(Dna)):
        for j in range(1, len(Dna)):
            if i + j >= len(Dna):
                temp = Dna[i:len(Dna)] + Dna[0: j - (len(Dna) - i)]
            else:
                temp = Dna[i:i+j]
            num = 0
            for r in range(0, len(temp)):
                num = num + mass[temp[r]]
            res.append(num)
    sum = 0
    for r in range(0, len(Dna)):
        sum = sum + mass[Dna[r]]
    res.append(sum)
    res.append(0)
    res.sort()
    return res

# mass of a peptide
def mass(peptide):
    infile = '/Users/yijia/Documents/current work/Bio/integer_mass_table.txt'
    mass = dict()
    with open(infile) as f:
        alph = f.readlines()
    Pair = [x.strip() for x in alph]
    for pair in Pair:
        mass[pair.split(' ')[0]] = int(pair.split(' ')[1])
    res = 0
    for i in range(0, len(peptide)):
        res = res + mass[peptide[i]]
    return res


#cycle for num
def cycle(num):
    res = list()
    for i in range(0, len(num)):
        for j in range(1, len(num)):
            if i + j >= len(num):
                temp = sum(num[i:len(num)]) + sum(num[0: j - (len(num) - i)])
            else:
                temp = sum(num[i:i + j])
            res.append(temp)
    if 0 not in res:
        res.append(0)
    res.append(sum(num))
    return sorted(res)

# whether a linear peptide is inconsistent with spectrum
def consistent(pep, spectrum):
    arr = [int(x) for x in pep.strip().split('-') if x != ""]
    linear = []
    for i in range(0, len(arr)):
        for j in range(1, len(arr) - i+1):
            temp = sum(arr[i: i+j])
            linear.append(temp)
    for m in linear:
        if m not in spectrum:
            return False
        if linear.count(m) > spectrum.count(m):
            return False
    return True

# cycle peptide sequence
def cyclopeptide(spectrum):
    #print consistent("103-137-71-131", spectrum)
    available = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    sortedspectrum = sorted(spectrum)
    maxValue = max(spectrum)
    #print maxValue
    res = set()
    peptide = []
    peptide.append("")
    while True:
        expend = []
        for p in peptide:
            #print p
            for s in spectrum:
                if s == 0 or s not in available:
                    continue
                if sum([int(x) for x in p.strip().split('-') if x != ""]) > maxValue:
                    continue
                if p != "":
                    temp = p + '-' + str(s)
                else:
                    temp = str(s)
                expend.append(temp)
        peptide = expend
        print len(expend)
        idx = len(peptide) - 1
        while idx >= 0:
            pp = peptide[idx]
            value = sum([int(x) for x in pp.strip().split('-') if x != ""])
            #print value
            if value == maxValue:
                if cycle([int(x) for x in pp.split('-')]) == sortedspectrum:
                    res.add(pp)
                peptide.remove(pp)
            elif value < maxValue:
                if not consistent(pp, spectrum):
                    peptide.remove(pp)
            elif value > maxValue:
                peptide.remove(pp)
            idx = idx - 1
        if len(peptide) == 0:
            break
    return sorted(res)

# 4.7 calculate score of peptide, spectrum
def score(peptide, spectrum):
    cycle = []
    for i in range(0, len(peptide)):
        for j in range(1, len(peptide)):
            if i + j >= len(peptide):
                temp = peptide[i:len(peptide)] + peptide[0: j - (len(peptide) - i)]
            else:
                temp = peptide[i:i + j]
            cycle.append(mass(temp))
    cycle.append(0)
    cycle.append(mass(peptide))
    cycle.sort()
    dif = 0
    set1 = set(cycle)
    set2 = set(spectrum)
    for n in set1:
        if n in spectrum:
            t = cycle.count(n) - spectrum.count(n)
        else:
            t = cycle.count(n)
        if t > 0:
            dif = dif + t
    return len(cycle) - dif

# 4.13.1 linear score
def linearScore(peptide, spectrum):
    cycle = []
    for i in range(0, len(peptide)):
        for j in range(1, len(peptide)-i+1):
            temp = peptide[i:i + j]
            cycle.append(mass(temp))
    cycle.append(0)
    #cycle.append(mass(peptide))
    cycle.sort()
    #print cycle
    dif = 0
    set1 = set(cycle)
    set2 = set(spectrum)
    for n in set1:
        if n in spectrum:
            t = cycle.count(n) - spectrum.count(n)
        else:
            t = cycle.count(n)
        if t > 0:
            dif = dif + t
    return len(cycle) - dif

# 4.13.3 trim
def trim(learerboard, spectrum, N):
    #print learerboard, N
    linear = []
    for j in range(0, len(learerboard)):
        peptide = learerboard[j]
        linear.append(linearScore(peptide, spectrum))
    learerboard = [x for (y, x) in sorted(zip(linear, learerboard))]
    learerboard.reverse()
    linear = sorted(linear, reverse=True)
    #print linear
    for j in range(N, len(learerboard)):
        if linear[j] < linear[N-1]:
            learerboard = learerboard[0:j]
            return learerboard
    return learerboard

# 4.7.10 linear score int
def scoreInt(pep, spectrum):
    peptide = [int(x) for x in pep.split("-") if x != ""]
    cycle = []
    for i in range(0, len(peptide)):
        for j in range(1, len(peptide)-i+1):
            temp = peptide[i:i + j]
            cycle.append(sum(temp))
    cycle.append(0)
    cycle.append(sum(peptide))
    cycle.sort()
    dif = 0
    set1 = set(cycle)
    set2 = set(spectrum)
    for n in set1:
        if n in spectrum:
            t = cycle.count(n) - spectrum.count(n)
        else:
            t = cycle.count(n)
        if t > 0:
            dif = dif + t
    res = len(cycle) - dif
    return res

#3print scoreInt('114-128-129-113',[0,99,113,114,128,227,257,299,355,356,370,371,484])

def scoreIntCycle(pep, spectrum):
    peptide = [int(x) for x in pep.split("-") if x != ""]
    cycle = []
    for i in range(0, len(peptide)):
        for j in range(1, len(peptide)):
            if i + j >= len(peptide):
                temp = peptide[i:len(peptide)] + peptide[0: j - (len(peptide) - i)]
            else:
                temp = peptide[i:i + j]
            cycle.append(sum(temp))
    cycle.append(0)
    cycle.append(sum(peptide))
    cycle.sort()
    dif = 0
    set1 = set(cycle)
    set2 = set(spectrum)
    for n in set1:
        if n in spectrum:
            t = cycle.count(n) - spectrum.count(n)
        else:
            t = cycle.count(n)
        if t > 0:
            dif = dif + t
    res = len(cycle) - dif
    return res

# trimInt for 4.7.8
def trimInt(learerboard, spectrum, N):
    #print learerboard, N
    linear = []
    for j in range(0, len(learerboard)):
        peptide = learerboard[j]
        linear.append(scoreInt(peptide, spectrum))
    learerboard = [x for (y, x) in sorted(zip(linear, learerboard))]
    learerboard.reverse()
    linear = sorted(linear, reverse=True)
    #print linear
    for j in range(N, len(learerboard)):
        if linear[j] < linear[N-1]:
            learerboard = learerboard[0:j]
            return learerboard
    #print "trim", learerboard
    return learerboard

#tryset = ['147-71','147-147', '147-113', '147-129', '147-147', '129-71', '129-113', '129-129', '129-147', '113-71', '113-113', '113-129', '113-147', '71-71', '71-113', '71-129', '71-147']
#print trimInt(tryset, [0,71,113,129,147,200,218,260,313,331,347,389,460], 5)

def changeFormat(sequence):
    infile = '/Users/yijia/Documents/current work/Bio/integer.txt'
    mass = dict()
    with open(infile) as f:
        alph = f.readlines()
    Pair = [x.strip() for x in alph]
    for pair in Pair:
        mass[pair.split(' ')[1]] = (pair.split(' ')[0])
    arr = sequence.strip().split('-')
    res = ""
    for x in arr:
        if x == "":
            continue
        res = res + mass[x]
    #print res
    return res

def changeback(sequence):
    infile = '/Users/yijia/Documents/current work/Bio/integer.txt'
    mass = dict()
    with open(infile) as f:
        alph = f.readlines()
    Pair = [x.strip() for x in alph]
    for pair in Pair:
        mass[pair.split(' ')[0]] = (pair.split(' ')[1])
    res = ""
    for i in range(0, len(sequence)):
        x = sequence[i]
        if res == "":
            res = res + mass[x]
        else:
            res = res + '-' + mass[x]
    #print res
    return res

def leadboard(spectrum, N):
    available = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    maxValue = max(spectrum)
    leader = set()
    leader.add("")
    best = ""
    while True:
        expend = set()
        for p in leader:
            for s in spectrum:
                if s == 0 or s not in available:
                    continue
                if sum([int(x) for x in p.strip().split('-') if x != ""]) > maxValue:
                    continue
                if p != "":
                    temp = p + '-' + str(s)
                else:
                    temp = str(s)
                expend.add(temp)
        leader = expend
        #print len(expend)
        idx = len(leader) - 1
        deduct = set()
        for pp in leader:
            deduct = set()
            value = sum([int(x) for x in pp.strip().split('-') if x != ""])
            if value == maxValue:
                if linearScore(changeFormat(pp), spectrum) > linearScore(changeFormat(best), spectrum):
                    #print linearScore(changeFormat(pp), spectrum), linearScore(changeFormat(best), spectrum)
                    best = pp
            elif value > maxValue:
                #leader.remove(pp)
                deduct.add(pp)
            idx = idx - 1
        leader = leader - deduct
        leadform = [changeFormat(x) for x in leader]
        leadform = trim(leadform, spectrum, N)
        leader = set([changeback(x) for x in leadform])
        if len(leader) == 0:
            break
    return best


'''
# leardboard cycleopeptide sequence 4.7.8
def leadSequence(spectrum, N):
    available = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    maxValue = max(spectrum)
    maxLen = maxValue/available[0]
    leadboard = set()
    leaderPeptide = ""
    leadboard.add("")
    while True:
        expend = []
        for p in leadboard:
            for s in spectrum:
                if s == 0 or s not in available:
                    continue
                if sum([int(x) for x in p.strip().split('-') if x != ""]) > maxValue:
                    continue
                if len(p.strip().split('-')) >= maxLen:
                    continue
                if p != "":
                    temp = p + '-' + str(s)
                else:
                    temp = str(s)
                expend.append(temp)
        print "expend", len(expend)
        leadboard = expend
        idx = len(leadboard) - 1
        while idx >= 0:
            peptide = leadboard[idx]
            value = sum([int(x) for x in peptide.strip().split('-') if x != ""])
            #print value
            if value == maxValue:
                #if scoreInt(peptide, spectrum) > scoreInt(leaderPeptide, spectrum):
                if linearScore(changeFormat(peptide), spectrum) > linearScore(changeFormat(peptide), spectrum):
                    leaderPeptide = peptide
            elif value > maxValue:
                leadboard.remove(peptide)
            idx = idx - 1
        #leadboard = trimInt(leadboard, spectrum, N)
        leadform = [changeFormat(x) for x in leadboard]
        leadform = trim(leadform, spectrum, N)
        leadboard = [changeback(x) for x in leadform]
        print len(leadboard)
        if len(leadboard) == 0:
            break
        #print leadboard
    return leaderPeptide
'''

# 4.9.4
def spectralConvolution(spectrum):
    res = []
    spectrum.sort()
    print spectrum
    for i in range(1, len(spectrum)):
        for j in range(0, i):
            res.append(spectrum[i] - spectrum[j])
    res.sort()
    idx = 0
    while res[idx] == 0:
        res.remove(0)
    return res

infile = '/Users/yijia/Documents/current work/Bio/dataset_38448_4.txt'
with open(infile) as f:
    Rna = f.readlines()
#peptide = Rna[0].strip().split(' ')
#N = int(Rna[0].strip())
#spectrum = [int(x) for x in Rna[1].split(" ")]
#spectrum = [0,71,113,129,147,200,218,260,313,331,347,389,460]
#N = 10
res = spectralConvolution([int(x) for x in Rna[0].split(" ")])
print res
outfile = open('/Users/yijia/Documents/current work/Bio/dataset_38448_4_res.txt', 'w')
for p in res:
    outfile.write(str(p))
    outfile.write(' ')
outfile.close()
