import sys
import random
import numpy as np
sys.setrecursionlimit(100000000) # 10000 is an example, try with different values

# 14.5.8
def HiddenPath(text, matrix):
    res = 0.5
    for i in range(1, len(text)):
        res = res*matrix[text[i]][text[i-1]]
    return res

# 14.5.10
def HMM(hidden, emission, matrix):
    res = 1
    for i in range(0, len(hidden)):
        res = res*matrix[emission[i]][hidden[i]]
    return res

# 14.6.7
def viterbi(text, transition, emission):
    res = ""
    states = np.array(["A", "B"])
    hidden = ["x", "y", "z"]
    place = dict()
    s = [[0 for j in range(0, len(text))] for i in range(0, len(states))]
    for i in range(0, len(states)):
        s[i][0] = 1.0/len(states)*emission[states[i]][text[0]]
    for i in range(1, len(text)):
        for j in range(0, len(states)):
            temp = []
            for k in range(0, len(states)):
                temp.append(s[k][i-1] * transition[states[k]][states[j]] * emission[states[j]][text[i]])
            s[j][i] = max(temp)
            place[(j, i)] = (temp.index(max(temp)), i-1)
    cur = ([m[-1] for m in s].index(max([m[-1] for m in s])), len(text) - 1)
    while cur[1] != 0:
        res = states[cur[0]] + res
        cur = place[cur]
    res = states[cur[0]] + res
    return res

# 14.7.3
def outcome(text, transition, emission):
    res = 0
    states = np.array(["A", "B", "C"])
    hidden = ["x", "y", "z"]
    place = dict()
    s = np.zeros((len(states), len(text)))
    for i in range(0, len(states)):
        s[i][0] = 1.0 / len(states) * emission[states[i]][text[0]]
    for i in range(1, len(text)):
        for j in range(0, len(states)):
            temp = []
            for k in range(0, len(states)):
                temp.append(s[k][i - 1] * transition[states[k]][states[j]] * emission[states[j]][text[i]])
            s[j][i] = sum(temp)
    return sum(s[:, len(text)-1])

h = "xyxzzxyxyy"
m = dict()
t = dict()
t["A"] = {"A": 0.641, "B": 359}
t["B"] = {"A": 0.729, "B": 0.271}
# t["C"] = {"A": 0.319, "B": 0.469, "C": 0.212}
m["A"] = {"x": 0.117, "y": 0.691, "z": 0.192}
m["B"] = {"x": 0.097, "y": 0.42, "z": 0.483}
# m["C"] = {"x": 0.67, "y": 0.314, "z": 0.016}
# print outcome(h, t, m)
print viterbi(h, t, m)