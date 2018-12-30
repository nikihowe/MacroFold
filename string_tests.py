import time
import csv
import math
import numpy as np

global base_dict
base_dict = {'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'U': 3, 'u': 3, 'T': 3, 't': 3}

global int_dict
int_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}

# Fill the 2x2 internal bulge tables (sequence : energy)
def loadInt222(zero_energies=False):
    global int222
    int222 = np.full((4,4,4,4,4,4,4,4), 1, dtype=np.float64)
    if zero_energies:
        filename = "zero_int22.csv"
    else:
        filename = "int22.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the first row
            if i == 0:
                continue
            else:
                int222[base_dict[row[0]],base_dict[row[1]],base_dict[row[2]],base_dict[row[3]],
                      base_dict[row[4]],base_dict[row[5]],base_dict[row[6]],base_dict[row[7]]] = math.exp(-1.6*float(row[8]))

# Fill the 2x2 internal bulge tables (sequence : energy)
def loadInt22(zero_energies=False):
    global int22
    int22 = {}
    if zero_energies:
        filename = "zero_int22.csv"
    else:
        filename = "int22.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the first row
            if i == 0:
                continue
            else:
                int22[str(row[0])+str(row[1])+str(row[2])+str(row[3])+
                      str(row[4])+str(row[5])+str(row[6])+str(row[7])] = math.exp(-1.6*float(row[8]))


endPenalty = {('AU'): 1, ('UA'): 1, ('GU'): 1, ('UG'): 1}
loopEndPenalty = {('AU'): 1, ('UA'): 1, ('GU'): 1, ('UG'): 1}
loopFirstMismatchEnergy = {('UU'): 1, ('GA'): 1, ('AG'): 1}
guClosureEnergy = {('GGGU'): 1}

def load():
    global s
    global t
    s = ""
    t = []
    for i in range(100000):
        s += int_dict[i%4]
        t.append(i%4)
    loadInt22()
    loadInt222()

def tuple_time():
    t0 = time.time()
    for i in range(1000000):
        a = ''.join(('a', 'b', 'c', 'd', str(i+3)))
        b = a[2]+a[3]
    t1 = time.time()
    print("tuple:", t1-t0)

def plus_time():
    t0 = time.time()
    for i in range(1000000):
        a = 'a'+'b'+'c'+'d'+str(i+3)
        b = a[2]+a[3]
    t1 = time.time()
    print("plus:", t1-t0)

def list_time():
    t0 = time.time()
    for i in range(1000000):
        a = ''.join(['a', 'b', 'c', 'd', str(i+3)])
        b = a[2]+a[3]
    t1 = time.time()
    print("list:", t1-t0)

def integer_rep():
    t0 = time.time()
    for i in range(1000000):
        a = [0, 1, 2, 3, 5]
        b = [a[2], a[3]]
    t1 = time.time()
    print("integer:", t1-t0)

def access_set():
    t0 = time.time()
    for i in range(99990):
        b = int22.get(s[i]+s[i+1]+s[i+2]+s[i+3]+s[i+4]+s[i+5]+s[i+6]+s[i+7], 1)
    t1 = time.time()
    print("set:", t1-t0)

def access_np():
    t0 = time.time()
    for i in range(99990):
        b = int222[t[i],t[i+1],t[i+2],t[i+3],t[i+4],t[i+5],t[i+6],t[i+7]]
    t1 = time.time()
    print("np:", t1-t0)

load()
access_set()
access_np()
