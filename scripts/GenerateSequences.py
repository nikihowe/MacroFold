# Nikolaus Howe August 2016

import numpy as np

def intToBase(n):
    intDict = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}
    return intDict.get(n, 'N')

def getComplement(base):
    compDict = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
    return compDict.get(base, 'N')

# Generate hairpins
def getHairpin(side_length, hairpin_length):
    sequence = ""
    # Fill in the hairpin
    for i in range(hairpin_length):
        sequence += intToBase(np.random.randint(0, 4))
    # Fill in to the right and left of the hairpin
    for i in range(side_length):
        nextbase = intToBase(np.random.randint(0, 4))
        sequence = nextbase + sequence + getComplement(nextbase)
    return sequence

# Generate internal loops and hairpins
def getLoop(hairpin_length, side2_length, top, bottom, side1_length):
    sequence = ""
    for i in range(hairpin_length):
        sequence += intToBase(np.random.randint(0, 4))
    # Fill in to the right and left of the hairpin
    for i in range(side2_length):
        nextbase = intToBase(np.random.randint(0, 4))
        sequence = nextbase + sequence + getComplement(nextbase)
    # Fill in internal loop
    for i in range(top):
        sequence = intToBase(np.random.randint(0, 4)) + sequence
    for i in range(bottom):
        sequence = sequence + intToBase(np.random.randint(0, 4))
    # Fill in last stacking parts
    for i in range(side1_length):
        nextbase = intToBase(np.random.randint(0, 4))
        sequence = nextbase + sequence + getComplement(nextbase)    
    return sequence

################
## Generators ##
################

# Generate 1000 different hairpin + internal loops as list
def generateLoops():
    output = []
    for i in range(300):
        output.append(getLoop(3, 3, 2, 1, 3))
    for i in range(300):
        output.append(getLoop(3, 3, 3, 0, 4))
    for i in range(300):
        output.append(getLoop(4, 5, 2, 2, 3))
    return output

# Generate 1000 different hairpins and save as list
def generateHairpins():
    output = []
    for i in range(300):
        output.append(getHairpin(1, 3))
    for i in range(300):
        output.append(getHairpin(3, 3))
    for i in range(300):
        output.append(getHairpin(4, 3))
    return output
