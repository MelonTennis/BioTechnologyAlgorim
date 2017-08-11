import numpy as np


# 15.1.15
def HMMProblem(theta, alpha, alignment):
    n = len(alignment[0])
    k = n * 3 + 3
    trans = [[0 for i in range(0, len(alignment) * n)] for j in range(0, k)]
    emission = [[0 for i in range(0, len(alpha))] for j in range(0, k)]
    states = []  # i, m, d
    for i in range(0, n):
        temp = [s[i] for s in alignment]
        frac = float(temp.count('-')) / k
        dele = []
        notdele = []
        for m in range(len(alignment)):
            if alignment[m][i] != '-':
                notdele.append(m)
            else:
                dele.append(m)
        if frac > theta:  # insert
            cur = [notdele, [], dele]
        else:
            cur = [[], notdele, dele]
        states.append(list(cur))
    idx = 0
    for s in len(states):
        if s != 0:
            prev_state = states[s - 1]
        else:
            prev_state = None
        if len(states[s][0]) != 0:  # insert
            if not prev_state or len(states[s - 1][0]) != 0:  # last is not insert
                # insert
                trans = calculate(prev_state, states[s], trans, idx, 'insert', s)
        else:  # match
            if not prev_state or len(states[s - 1][0]) == 0:  # last is not insert
                trans = calculate(prev_state, states[s], trans, idx, 'match', s)
            else:
                trans = calculate(prev_state, states[s], trans, idx, 'transform', s)
            idx += 1

    idx = 0
    for s in len(states):
        if len(states[s][0]) != 0:  # insert
            emission = emit(emission, states[s], idx, alignment, 'insert', alpha)
        else:
            emission = emit(emission, states[s], idx, alignment, 'match', alpha)
            idx += 1
    return trans, emission


def calculate(prev_state, state, trans, idx, form, pos):
    if form == 'insert':
        if not prev_state:
            trans[0][1] = float(len(state[0])) / (len(state[0]) + len(state[1]) + len(state[2]))
        else:
            # match
            if len(prev_state[1]) > 0:
                trans[1 + 3 * (idx - 1) + 1][1 + 3 * (idx - 1) + 3] = float(
                    len([x for x in state[0] if x in prev_state[1]])) / len(prev_state[1])
            # delete
            if len(prev_state[2]) > 0:
                trans[1 + 3 * (idx - 1) + 2][1 + 3 * (idx - 1) + 3] = float(
                    len([x for x in state[0] if x in prev_state[2]])) / len(prev_state[2])
    elif form == 'match':
        if not prev_state:
            trans[0][2] = float(len(state[1])) / (len(state[0]) + len(state[1]) + len(state[2]))
        else:
            # match
            if len(prev_state[1]) > 0:
                trans[1 + 3 * (idx - 1) + 1][1 + 3 * (idx - 1) + 1] = float(
                    len([x for x in state[1] if x in prev_state[1]])) / len(prev_state[1])
                trans[1 + 3 * (idx - 1) + 1][1 + 3 * (idx - 1) + 2] = float(
                    len([x for x in state[2] if x in prev_state[1]])) / len(prev_state[1])
            # delete
            if len(prev_state[2]) > 0:
                trans[1 + 3 * (idx - 1) + 2][1 + 3 * (idx - 1) + 1] = float(
                    len([x for x in state[1] if x in prev_state[2]])) / len(prev_state[2])
                trans[1 + 3 * (idx - 1) + 2][1 + 3 * (idx - 1) + 2] = float(
                    len([x for x in state[2] if x in prev_state[2]])) / len(prev_state[2])
    else:
        # insert
        if len(prev_state[0]) > 0:
            trans[1 + 3 * (idx - 1)][1 + 3 * (idx - 1) + 1] = float(
                len([x for x in state[1] if x in prev_state[0]])) / len(prev_state[0])
            trans[1 + 3 * (idx - 1)][1 + 3 * (idx - 1) + 2] = float(
                len([x for x in state[2] if x in prev_state[0]])) / len(prev_state[0])
    return trans

def emit(emission, state, idx, alignment, form, alpha):
    if form == 'insert':
        for i in range(0, len(alpha)):
            emission[idx][i] = float([x[idx] for x in alignment].count(state[i])) / len(state[0])
    else:
        for i in range(0, len(alpha)):
            emission[idx][i] = float([x[idx] for x in alignment].count(state[i])) / len(state[1])

def HMMProfile(theta, alphabet, alignmnet):
    graph, diction, num = prePorcess(theta, alphabet, alignmnet)
    trans, states = countTrans(graph, num)
    emit = countEmit(diction, alphabet, states)
    return trans, emit, states

def prePorcess(theta, alphabet, alignment):
    cnt = 0
    num = len(alignment[0])
    emits = dict()
    k = len(alignment[0])
    graph = [[] for i in range(0, len(alignment))]
    for i in range(0, k):
        curAlign = [x[i] for x in alignment]
        frac = float(curAlign.count('-'))/len(alignment)
        if frac > theta:
            # insert
            num -= 1
            for j in range(0, len(alignment)):
                if alignment[j][i] != "-":
                    cur_state = "I" + str(cnt)
                    graph[j].append(cur_state)
                    if cur_state not in emits.keys():
                        emits[cur_state] = []
                    emits[cur_state].append(alignment[j][i])
        else:
            cnt += 1
            for j in range(0, len(alignment)):
                if alignment[j][i] != "-":
                    cur_state = "M" + str(cnt)
                else:
                    cur_state = "D" + str(cnt)
                graph[j].append(cur_state)
                if cur_state not in emits.keys():
                    emits[cur_state] = []
                emits[cur_state].append(alignment[j][i])
    # print graph
    return graph, emits, num

def countTrans(graph, num):
    states = ["S", "I0"]
    function = dict()
    function["M"] = lambda x: 3*x - 1
    function["D"] = lambda x: 3*x
    function["I"] = lambda x: 3*x + 1
    matrix = [[0 for i in range(0, num*3 + 3)] for j in range(0, num*3 + 3)]
    for i in range(0, len(graph)):
        for j in range(0, len(graph[i])):
            state = graph[i][j][0]
            idx = int(graph[i][j][1])
            if j == 0:
                if state == "I":
                    matrix[0][1] += 1
                elif state == "M":
                    matrix[0][2] += 1
                elif state == "D":
                    matrix[0][3] += 1
            else:
                pre_state = graph[i][j-1][0]
                pre_idx = int(graph[i][j-1][1])
                # print function[pre_state](pre_idx), function[state](idx)
                # print matrix
                matrix[function[pre_state](pre_idx)][function[state](idx)] += 1
            if j == len(graph[i]) - 1:
                matrix[function[state](idx)][-1] += 1
    for i in range(1, num + 1):
        states.append("M" + str(i))
        states.append("D" + str(i))
        states.append("I" + str(i))
    states.append("E")
    print states
    assert len(states) == len(matrix)
    for l in range(0, len(matrix)):
        hold = sum(matrix[l])
        if hold != 0:
            matrix[l] = [(float(x)/hold) if x != 0 else x for x in matrix[l]]
    # for j in range(0, len(matrix)):
    #     matrix[j] = [states[j]] + matrix[j]
    # matrix = [states] + matrix
    print states
    return matrix, states

def countEmit(diction, alphabet, states):
    function = dict()
    function["M"] = lambda x: 3 * x - 1
    function["D"] = lambda x: 3 * x
    function["I"] = lambda x: 3 * x + 1
    matrix = [[0 for i in range(0, len(alphabet))] for j in range(0, len(states))]
    for key in diction.keys():
        state = key[0]
        idx = int(key[1])
        row = function[state](idx)
        for a in alphabet:
            matrix[row][alphabet.index(a)] = diction[key].count(a)
    for l in range(0, len(matrix)):
        hold = sum(matrix[l])
        if hold != 0:
            matrix[l] = [(float(x)/hold) if x != 0 else x for x in matrix[l]]
    # for j in range(0, len(matrix)):
    #     matrix[j] = [states[j]] + matrix[j]
    # matrix = [alphabet] + matrix
    return matrix

# 15.2.5
def pseuHMM(theta, rate, alpha, alignment):
    function = dict()
    function["M"] = lambda x: 3 * x - 1
    function["D"] = lambda x: 3 * x
    function["I"] = lambda x: 3 * x + 1
    trans, emit, states = HMMProfile(theta, alpha, alignment)
    for i in range(0, len(trans)):
        idx = (i + 1)/3
        for j in range((idx + 1)*3-2, min((idx + 1)*3 + 1, len(trans[i]))):
            if trans[i][j] == 0:
                trans[i][j] = rate
    for i in range(1, len(emit) - 1):
        if i%3 == 0:
            continue
        else:
            emit[i] = [x if x != 0 else rate for x in emit[i]]
    for l in range(0, len(trans)):
        hold = sum(trans[l])
        if hold != 0:
            trans[l] = ["%.3f"%(float(x)/hold) if x != 0 else x for x in trans[l]]
    for j in range(0, len(trans)):
        trans[j] = [states[j]] + trans[j]
    trans = [states] + trans
    for l in range(0, len(emit)):
        hold = sum(emit[l])
        if hold != 0:
            emit[l] = ["%.3f"%(float(x)/hold) if x != 0 else x for x in emit[l]]
    for j in range(0, len(emit)):
        emit[j] = [states[j]] + emit[j]
    emit = [alpha] + emit
    return trans, emit

# 15.4.4
def HMMParameter(x, alpha, path, states):
    trans = [[0 for i in range(0, len(states))] for j in range(0, len(states))]
    emiss = [[0 for i in range(0, len(alpha))] for j in range(0, len(states))]
    countStates = [path[0:-1].count(m) for m in states]
    for i in range(0, len(path) - 1):
        trans[states.index(path[i])][states.index(path[i + 1])] += 1
    for i in range(0, len(states)):
        if sum(trans[i]) != 0:
            trans[i] = [float(sl) / countStates[i] for sl in trans[i]]
        else:
            trans[i] = [1.0 / len(states) for k in range(0, len(states))]
            # trans[i] = [states[i]] + trans[i]
    countEmit = [path.count(m) for m in states]
    for i in range(0, len(x)):
        emiss[states.index(path[i])][alpha.index(x[i])] += 1
    for i in range(0, len(states)):
        if sum(emiss[i]) != 0:
            emiss[i] = [(float(sl) / countEmit[i]) for sl in emiss[i]]
        else:
            emiss[i] = [1.0 / len(alpha) for k in range(0, len(alpha))]
            # emiss[i] = [states[i]] + emiss[i]
    return trans, emiss


# 15.4.8
def viterbiIter(iters, x, alpha, states, trans, emission):
    for iter in range(0, iters):
        path = viterbi(x, trans, emission, states, alpha)
        trans, emission = HMMParameter(x, alpha, path, states)
    for i in range(0, len(states)):
        trans[i] = [states[i]] + ["%.3f" % k for k in trans[i]]
        emission[i] = [states[i]] + ["%.3f" % k for k in emission[i]]
    trans = [states] + trans
    emission = [alpha] + emission
    return trans, emission


def viterbi(text, transition, emission, states, alpha):
    res = ""
    states = np.array(states)
    place = dict()
    s = [[0 for j in range(0, len(text))] for i in range(0, len(states))]
    for i in range(0, len(states)):
        s[i][0] = 1.0 / len(states) * emission[i][alpha.index(text[0])]
    for i in range(1, len(text)):
        for j in range(0, len(states)):
            temp = []
            for k in range(0, len(states)):
                temp.append(s[k][i - 1] * transition[k][j] * emission[j][alpha.index(text[i])])
            s[j][i] = max(temp)
            place[(j, i)] = (temp.index(max(temp)), i - 1)
    cur = (np.argmax([m[-1] for m in s]), len(text) - 1)
    while cur[1] != 0:
        res = states[cur[0]] + res
        cur = place[cur]
    res = states[cur[0]] + res
    return res

# 15.5.5
def forwardBackward(x, alpha, states, trans, emit):
    forward = [[0 for j in range(len(states))] for i in range(0, len(x))]
    backward = [[0 for j in range(len(states))] for i in range(0, len(x))]
    weight = [[0 for j in range(0, len(states))] for i in range(len(x))]
    prob = [1.0 / len(states) for i in range(0, len(states))]
    for i in range(0, len(states)):
        forward[0][i] = prob[i] * emit[i][alpha.index(x[0])]
        backward[-1][i] = 1
    for i in range(1, len(x)):
        for j in range(0, len(states)):
            for k in range(0, len(states)):
                forward[i][j] += forward[i - 1][k] * trans[k][j]
            forward[i][j] *= emit[j][alpha.index(x[i])]
    sink = sum([w[-1] for w in forward])
    print sink
    # print forward
    for i in range(len(x) - 2, -1, -1):
        for j in range(0, len(states)):
            for k in range(0, len(states)):
                backward[i][j] += trans[j][k] * emit[k][alpha.index(x[i + 1])] * backward[i + 1][k]
    print backward
    for i in range(0, len(x)):
        for j in range(0, len(states)):
            weight[i][j] = "%.3f" % (forward[i][j] * backward[i][j] / float(sink))
    weight = [states] + weight
    return weight

# 15.6.5
def forwardBackward1(x, alpha, states, trans, emit):
    forward = [[0 for j in range(len(states))] for i in range(0, len(x))]
    backward = [[0 for j in range(len(states))] for i in range(0, len(x))]
    weight = [[0 for j in range(0, len(states))] for i in range(len(x))]
    prob = [1.0 / len(states) for i in range(0, len(states))]
    para = [[[0 for i in range(0, len(x) - 1)] for j in range(0, len(states))] for k in range(0, len(states))]
    for i in range(0, len(states)):
        forward[0][i] = prob[i] * emit[i][alpha.index(x[0])]
        backward[-1][i] = 1
    for i in range(1, len(x)):
        for j in range(0, len(states)):
            for k in range(0, len(states)):
                forward[i][j] += forward[i - 1][k] * trans[k][j]
            forward[i][j] *= emit[j][alpha.index(x[i])]
    sink = sum(forward[-1])
    for i in range(len(x) - 2, -1, -1):
        for j in range(0, len(states)):
            for k in range(0, len(states)):
                backward[i][j] += trans[j][k] * emit[k][alpha.index(x[i + 1])] * backward[i + 1][k]
    for i in range(0, len(x)):
        for j in range(0, len(states)):
            weight[i][j] = forward[i][j] * backward[i][j] / float(sink)
    for i in range(0, len(states)):
        for j in range(0, len(states)):
            para[i][j] = []
            for k in range(0, len(x)-1):
                w = trans[i][j] * emit[j][alpha.index(x[k+1])]
                para[i][j].append(forward[k][i] * backward[k+1][j] * w)
    return para, weight

def BWLearning(iter, x, alpha, states, trans, emit):
    for i in range(0, iter):
        para, weight = forwardBackward1(x, alpha, states, trans, emit)
        trans = update(para, weight, trans, "t", states, alpha, x)
        emit = update(para, weight, emit, "e", states, alpha, x)
    return trans, emit

def update(para, weight, ori, form, states, alpha, x):
    if form == 't':
        for i in range(0, len(states)):
            for j in range(0, len(states)):
                ori[i][j] = sum(para[i][j])
        for i in range(0, len(states)):
            ori[i] = [float(kk)/sum(ori[i]) for kk in ori[i]]
    if form == 'e':
        renew = [[0 for i in range(0, len(alpha))] for j in range(0, len(states))]
        for i in range(0, len(x)):
            for j in range(0, len(states)):
                renew[j][alpha.index(x[i])] += weight[i][j]
        for i in range(0, len(states)):
            renew[i] = [float(kk)/sum(renew[i]) for kk in renew[i]]
        ori = renew
    return ori

# iter = 100
# x = "zyxxxxyxzz"
# alpha = ["x", "y", "z"]
# state = ["A", "B"]
# trans = [[0.911, 0.089], [0.228, 0.772]]
# emit = [[0.356, 0.191, 0.453], [0.040, 0.467, 493]]
# forwardBackward(x, alpha, state, trans, emit)
# tran, emit = BWLearning(iter, x, alpha, state, trans, emit)
# for i in range(0, len(state)):
#     tran[i] = [state[i]] + ["%.3f"% r for r in tran[i]]
# tran = [state] + tran
# for i in range(0, len(state)):
#     emit[i] = [state[i]] + ["%.3f"%r for r in emit[i]]
# emit = [alpha] + emit
# print tran
# print emit
# weight = forwardBackward(x, alpha, state, trans, emit)
# print weight

infile = '/Users/yijia/Documents/current_work/Bio/file.txt'
with open(infile) as f:
    input = f.readlines()
theta = float(input[0].strip().split(" ")[0])
rate = float(input[0].strip().split(" ")[1])
alpha = input[2].strip().split(" ")
alignments = [x.strip() for x in input[4:]]
trans, emit = pseuHMM(theta, rate, alpha, alignments)
print trans
print emit
outfile = open('/Users/yijia/Documents/current_work/Bio/file_res.txt', 'w')
for tr in trans:
    for r in tr:
        outfile.write(str(r) + "\t")
    outfile.write("\n")
outfile.write("--------")
outfile.write("\n")
for er in emit:
    for r in er:
        outfile.write(str(r) + "\t")
    outfile.write("\n")
