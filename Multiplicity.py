# MacroFold.py
# (c) 2016-2017 Nikolaus Howe
# by Niki Howe,    
#                  Summer/Fall/Winter
#               (July->January) 2016/2017

# Imports
import numpy as np
import math
import sys
import csv
import random

# For testing with UNAFold
import subprocess

# The temperature (in Kelvin)
global T
T = 273.15 + 37

# The gas constant (kcal/mol*K) ## R=k_B*N_A=(1.38e-23 J/K)/(cal/4.184J)(6.022e23)  ##
global r
r = 0.0019872036

# The beta constant
global beta
beta = 1./(r*T)

# get_all_structures
def get_structs(sequence, only_num=False, write=True):

    # Make sure it is upper-case
    sequence = sequence.upper()

    # List of allowed pairs
    # Assumes only capital letters and no T's
    legal_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}

    # Make a new file which contains the sequence
    with open("./sequence.seq", 'w') as myfile:
        myfile.write(sequence+'\n')

    # Get all the folds of that sequence using crumple
    all_folds = subprocess.getoutput("./crumple/crumple -i sequence.seq")

    # Put them in a list and get rid of the last two lines of text
    all_folds = set(all_folds.split()[:-2])

    # The two places where dangles
    # can appear are in the external 
    # loop and in multibranches
    
    # External loop
    old_folds = set() 
    while all_folds != old_folds:
        old_folds = set(all_folds)
        for fold in old_folds:
            i = 0
            while i < len(fold) - 1:
                if fold[i] == '(':
                    if i > 0 and fold[i-1] == '.':
                        new_fold = fold[:i-1] + '<' + fold[i:]
                        all_folds.add(new_fold)
                    i = get_match(i, fold)
                if fold[i] == ')':
                    if i < len(fold)-1 and fold[i+1] == '.':
                        new_fold = fold[:i+1] + '>' + fold[i+2:]
                        all_folds.add(new_fold)
                i += 1
    
    # Multibranch loops
    old_folds = set()
    while all_folds != old_folds:
        old_folds = set(all_folds)
        for fold in old_folds:
            for i, char in enumerate(fold):
                if char == '(':
                    num = get_num_structures(i, fold)
                    if num >= 2:
                        j = get_match(i, fold)
                        #print(fold, i, j)
                        #if fold[i+1] == '.':
                            #new_fold = fold[:i+1] + '>' + fold[i+2:]
                            #all_folds.add(new_fold)
                        #if fold[j-1] == '.':
                            #new_fold = fold[:j-1] + '<' + fold[j:]
                        k = i+1
                        #print(i, k, j)
                        while k < j:
                            if fold[k] == '.':
                                if fold[k-1] == '(' or fold[k-1] == ')':
                                    new_fold = fold[:k] + '>' + fold[k+1:]
                                    all_folds.add(new_fold)
                                if fold[k+1] == '(' or fold[k+1] == ')':
                                    new_fold = fold[:k] + '<' + fold[k+1:]
                                    all_folds.add(new_fold)
                            if fold[k] == '(':
                                k = get_match(k, fold)
                            k += 1

    # Remove the folds with a tstack which is not allowed
    # (namely, if the stacking bases can pair with each other)
    old_folds = set(all_folds)
    for fold in old_folds:
        for i, char in enumerate(fold):
            if char == '(':
                if i > 0 and fold[i-1] == '<':
                    j = get_match(i, fold)
                    if j < len(fold) - 1 and fold[j+1] == '>':
                        if sequence[i-1]+sequence[j+1] in legal_pairs:
                            if fold in all_folds: # because we might have removed it in the meantime
                                all_folds.remove(fold)
                
                if fold[i+1] == '>':
                    j = get_match(i, fold)
                    if fold[j] == ')':
                        if fold[j-1] == '<':
                            if sequence[i+1]+sequence[j-1] in legal_pairs:
                                if fold in all_folds:
                                    all_folds.remove(fold)

    #for fold in old_folds - all_folds:
        #print(fold)
    #for fold in all_folds:
        #print(fold)

    # For every fold, make a pair_dict (right AND left-looking)
    # and a dangle_list, which we will use to make
    # the extended .ct files
    all_folds = sorted(list(all_folds))

    if only_num:
        return len(all_folds)
    
    structure_list = []
    
    for fold in all_folds:
        pair_dict = {}
        dangle_list = []
        for i, char in enumerate(fold):
            if char == '(':
                j = get_match(i, fold)
                pair_dict[i] = j
                pair_dict[j] = i
            if char == '<':
                dangle_list.append(i)
            if char == '>':
                dangle_list.append(i-1)
        structure_list.append((pair_dict, dangle_list))
    #for element in structure_list:
        #print(element)

    if write:
        # Write these folds to an extended .ct file
        for i, fold in enumerate(all_folds):
            #if (i+1)%2000 == 0:
                #print(i+1, fold)
            with open("./folds/fold{}.txt".format(i+1), "w") as text_file:
                # Get the corresponding structure
                pairs = structure_list[i][0]
                dangles = structure_list[i][1]
                
                #print(pairs, dangles)

                # Make lists for all the stacking
                left_stacking = []
                right_stacking = []

                # Stacking from pairing
                for pos, char in enumerate(fold):
                    if char == '(':
                        j = get_match(pos, fold)
                        num = get_num_structures(pos, fold)
                        if num == 1:
                            #print("one", fold[pos: j+1])
                            # normal stacking
                            if fold[pos+1] == '(' and fold[j-1] == ')':
                                left_stacking.append(pos+1)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-1)
                            # 0x1 bulge loop
                            elif fold[pos+1] == '(' and fold[j-1] == '.' and fold[j-2] == ')':
                                left_stacking.append(pos+1)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-2)
                            # 1x0 bulge loop
                            elif fold[pos+1] == '.' and fold[pos+2] == '(' and fold[j-1] == ')':
                                left_stacking.append(pos+2)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-1)

                # Stacking from dangles
                for d in dangles:
                    left_stacking.append(d+1)
                    right_stacking.append(d)

                #print("left stacking ", left_stacking)
                #print("right stacking", right_stacking)

                # Write the first line of the .ct file
                text_file.write(str(len(sequence))+'\t'+"fold{}".format(i+1)+'\n')

                # Write the other lines
                for line, char in enumerate(sequence):
                    # The first four columns
                    text_file.write(str(line+1)+'\t'+char+'\t'+str(line)+'\t'+str((line+2)%(len(sequence)+1))+'\t')

                    # The pairing column
                    if fold[line] == '(' or fold[line] == ')':
                        text_file.write(str(pairs[line]+1)+'\t')
                    else:
                        text_file.write(str(0)+'\t')

                    # Same as first column
                    text_file.write(str(line+1)+'\t')

                    # Stack-to-the-left column
                    if line > 0 and line in left_stacking:
                        if fold[line-1] == '(' or fold[line-1] == ')' or fold[line-1] == '<' or fold[line-1] == '>':
                            text_file.write(str(line)+'\t')
                        elif fold[line-1] == '.':
                            text_file.write(str(line-1)+'\t')
                        else:
                            raise Exception("This is not supposed to happen")
                    else:
                        text_file.write(str(0)+'\t')
                    
                    # Stack-to-the-right column
                    if line < len(fold) - 1 and line in right_stacking:
                        #print(line)
                        if fold[line+1] == ')' or fold[line+1] == '(' or fold[line+1] == '<' or fold[line+1] == '>':
                            text_file.write(str(line+2)+'\n')
                        elif fold[line+1] == '.':
                            text_file.write(str(line+3)+'\n')
                        else:
                            raise Exception("This is not supposed to happen")
                    else:
                        text_file.write(str(0)+'\n')
    # Returns the number of folds
    return(len(all_folds), all_folds, structure_list)


# Given a place to start looking, returns the number of structures inside this one
# Can be arbitrarily large
def get_num_structures(index, seq):

    # Get the closing index for this structure
    closing_index = get_match(index, seq) 
    #if index == 4 and closing_index == 20: print(seq, index, closing_index)
    next_index = index + 1
    num_structures = 0
    while next_index < closing_index:
        #if index == 4 and closing_index == 20: print(next_index)
        #if index == 4 and closing_index == 20: print(seq[next_index])
        if seq[next_index] == '(':
            num_structures += 1
            next_index = get_match(next_index, seq) # used to read get_match(next_index, seq) + 1. Fixed on Jan 16.
            #print("num_structures", num_structures)
        next_index += 1
    return num_structures


# index gives the index of the left paren we are looking
# for a match of; string gives the whole string in 
# which we are searching
def get_match(index, string, left_bracket='(', right_bracket=')'):

    # We keep track of how many parens we still
    # need to match
    count = 0

    # Add on left parens, subtract right ones
    # When there are zero left to match,
    # we have matched the original
    for i, char in enumerate(string[index:]):
        if char == left_bracket:
            count += 1
        if char == right_bracket:
            count -= 1
        if count <= 0:
            return index+i

#print(get_num_structures(0, "((...)(...).)"))

# Check the number of folds in a randomly generated sequence of length n
def generate_sequence(n):
    # The list of bases to sample from
    bases = ['A', 'C', 'G', 'U']
    # Sample with replacement with uniform weights
    result = ""
    for i in range(n):
        result += random.choice(bases)
    return result

# Given a sequence and a fold (in dot/bracket notation)
# save the .ct file
def make_ct(sequence, fold):

    # For every fold, make a pair_dict (right AND left-looking)
    # and a dangle_list, which we will use to make
    # the extended .ct files
    all_folds = sorted(list(all_folds))

    if only_num:
        return len(all_folds)
    
    structure_list = []
    
    for fold in all_folds:
        pair_dict = {}
        dangle_list = []
        for i, char in enumerate(fold):
            if char == '(':
                j = get_match(i, fold)
                pair_dict[i] = j
                pair_dict[j] = i
            if char == '<':
                dangle_list.append(i)
            if char == '>':
                dangle_list.append(i-1)
        structure_list.append((pair_dict, dangle_list))
    #for element in structure_list:
        #print(element)

    if write:
        # Write these folds to an extended .ct file
        for i, fold in enumerate(all_folds):
            #if (i+1)%2000 == 0:
                #print(i+1, fold)
            with open("./folds/fold{}.txt".format(i+1), "w") as text_file:
                # Get the corresponding structure
                pairs = structure_list[i][0]
                dangles = structure_list[i][1]
                
                #print(pairs, dangles)

                # Make lists for all the stacking
                left_stacking = []
                right_stacking = []

                # Stacking from pairing
                for pos, char in enumerate(fold):
                    if char == '(':
                        j = get_match(pos, fold)
                        num = get_num_structures(pos, fold)
                        if num == 1:
                            #print("one", fold[pos: j+1])
                            # normal stacking
                            if fold[pos+1] == '(' and fold[j-1] == ')':
                                left_stacking.append(pos+1)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-1)
                            # 0x1 bulge loop
                            elif fold[pos+1] == '(' and fold[j-1] == '.' and fold[j-2] == ')':
                                left_stacking.append(pos+1)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-2)
                            # 1x0 bulge loop
                            elif fold[pos+1] == '.' and fold[pos+2] == '(' and fold[j-1] == ')':
                                left_stacking.append(pos+2)
                                left_stacking.append(j)
                                right_stacking.append(pos)
                                right_stacking.append(j-1)

                # Stacking from dangles
                for d in dangles:
                    left_stacking.append(d+1)
                    right_stacking.append(d)

                #print("left stacking ", left_stacking)
                #print("right stacking", right_stacking)

                # Write the first line of the .ct file
                text_file.write(str(len(sequence))+'\t'+"fold{}".format(i+1)+'\n')

                # Write the other lines
                for line, char in enumerate(sequence):
                    # The first four columns
                    text_file.write(str(line+1)+'\t'+char+'\t'+str(line)+'\t'+str((line+2)%(len(sequence)+1))+'\t')

                    # The pairing column
                    if fold[line] == '(' or fold[line] == ')':
                        text_file.write(str(pairs[line]+1)+'\t')
                    else:
                        text_file.write(str(0)+'\t')

                    # Same as first column
                    text_file.write(str(line+1)+'\t')

                    # Stack-to-the-left column
                    if line > 0 and line in left_stacking:
                        if fold[line-1] == '(' or fold[line-1] == ')' or fold[line-1] == '<' or fold[line-1] == '>':
                            text_file.write(str(line)+'\t')
                        elif fold[line-1] == '.':
                            text_file.write(str(line-1)+'\t')
                        else:
                            raise Exception("This is not supposed to happen")
                    else:
                        text_file.write(str(0)+'\t')
                    
                    # Stack-to-the-right column
                    if line < len(fold) - 1 and line in right_stacking:
                        #print(line)
                        if fold[line+1] == ')' or fold[line+1] == '(' or fold[line+1] == '<' or fold[line+1] == '>':
                            text_file.write(str(line+2)+'\n')
                        elif fold[line+1] == '.':
                            text_file.write(str(line+3)+'\n')
                        else:
                            raise Exception("This is not supposed to happen")
                    else:
                        text_file.write(str(0)+'\n')
    # Returns the number of folds
    return(len(all_folds), all_folds, structure_list)
# Check if a variable has been assigned a value yet
def exists(x):
    return (x in locals() or x in globals())

# Fill all the energy tables
def fillEnergyTables(no_dangles=False, zero_energies=False, force=False):
    if force:
        loadTerminalMismatchHEnergy(no_dangles, zero_energies)
        loadTerminalMismatchEEnergy(no_dangles, zero_energies)
        loadTerminalMismatchMEnergy(no_dangles, zero_energies)
        loadDangle3Energy(no_dangles, zero_energies)
        loadDangle5Energy(no_dangles, zero_energies)
        loadBulgeLoopEnergy(zero_energies)
        loadInternalLoopEnergy(zero_energies)
        loadStackEnergy(zero_energies)
        loadHairpinEnergy(zero_energies)
        loadSpecialLoops(zero_energies)
        loadCLoops(zero_energies)
        loadAsymmetry(zero_energies)
        loadInt11(zero_energies)
        loadInt21(zero_energies)
        loadInt22(zero_energies)
        load_penalty_dicts(zero_energies)
        load_multibranch_penalties(zero_energies)
    else:
        if not exists("terminalMismatchHEnergy"): loadTerminalMismatchHEnergy(no_dangles, zero_energies)
        if not exists("terminalMismatchEEnergy"): loadTerminalMismatchEEnergy(no_dangles, zero_energies)
        if not exists("terminalMismatchMEnergy"): loadTerminalMismatchMEnergy(no_dangles, zero_energies)
        if not exists("dangle3Energy"): loadDangle3Energy(no_dangles, zero_energies)
        if not exists("dangle5Energy"): loadDangle5Energy(no_dangles, zero_energies)
        if not exists("bulgeLoopEnergy"): loadBulgeLoopEnergy(zero_energies)
        if not exists("internalLoopEnergy"): loadInternalLoopEnergy(zero_energies)
        if not exists("stackEnergy"): loadStackEnergy(zero_energies)
        if not exists("hairpinEnergy"): loadHairpinEnergy(zero_energies)
        if not exists("specialLoops"): loadSpecialLoops(zero_energies)
        if not exists("CLoops"): loadCLoops(zero_energies)
        if not exists("asymmetry"): loadAsymmetry(zero_energies)
        if not exists("int11"): loadInt11(zero_energies)
        if not exists("int21"): loadInt21(zero_energies)
        if not exists("int22"): loadInt22(zero_energies)
        if not exists("endPenalty"): load_penalty_dicts(zero_energies)
        if not exists("ea"): load_multibranch_penalties(zero_energies)

# Fill in the dangle3 table
# Note that the table is organized as: 5PrimeBase, 3PrimeBase, DangleBase, Energy
# where the 5PrimeBase (i) and 3PrimeBase (j) are paired, with DangleBase dangling off
# to the left of the 3PrimeBase (at j-1)
def loadDangle3Energy(no_dangle=False, zero_energies=False):
    global dangle3Energy
    dangle3Energy = {}
    if no_dangle:
        filename = "inf_dangle3.csv"
    elif zero_energies:
        filename = "zero_dangle3.csv"
    else:
        filename = "dangle3.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                dangle3Energy[str(row[0])+str(row[1])+str(row[2])] = math.exp(-beta*float(row[3]))
                
# Fill in the dangle5 table
# Note that the table is organized as: 5PrimeBase, 3PrimeBase, DangleBase, Energy
# where the 5PrimeBase (i) and 3PrimeBase (j) are paired, with DangleBase dangling off
# to the right of the 5PrimeBase (at i+1)
def loadDangle5Energy(no_dangle=False, zero_energies=False):
    global dangle5Energy
    dangle5Energy = {}
    if no_dangle:
        filename = "inf_dangle5.csv"
    elif zero_energies:
        filename = "zero_dangle5.csv"
    else:
        filename = "dangle5.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                dangle5Energy[str(row[0])+str(row[1])+str(row[2])] = math.exp(-beta*float(row[3]))
    
# Fill in the bulge loop energy dict (length : Energy)
def loadBulgeLoopEnergy(zero_energies=False):
    global bulgeLoopEnergy
    bulgeLoopEnergy = {}
    if zero_energies:
        filename = "zero_bulge.csv"
    else:
        filename = "bulge.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                bulgeLoopEnergy[int(row[0])] = math.exp(-beta*float(row[1]))

# Fill in the internal loop energy dict (length : Energy)
def loadInternalLoopEnergy(zero_energies=False):
    global internalLoopEnergy
    internalLoopEnergy = {}
    if zero_energies:
        filename = "zero_internal.csv"
    else:
        filename = "internal.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                internalLoopEnergy[int(row[0])] = math.exp(-beta*float(row[1]))

# Fill in the stack energy dict (5primeOuter5primeInner3primeInner3primeOuter : Energy)
def loadStackEnergy(zero_energies=False):
    global stackEnergy
    stackEnergy = {}
    if zero_energies:
        filename = "zero_stack.csv"
    else:
        filename = "stack.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                stackEnergy[str(row[0])+str(row[1])+str(row[2])+str(row[3])] = math.exp(-beta*float(row[4]))
                
# Fill in the stack terminal mismatch energy dict (5primeOuter5primeInner3primeInner3primeOuter : Energy)
def loadTerminalMismatchHEnergy(no_dangle=False, zero_energies=False):
    global terminalMismatchHEnergy
    terminalMismatchHEnergy = {}
    if no_dangle:
        filename = "inf_tstackh.csv"
    elif zero_energies:
        filename = "zero_tstackh.csv"
    else:
        filename = "tstackh.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                terminalMismatchHEnergy[str(row[0])+str(row[1])+str(row[2])+str(row[3])] = math.exp(-beta*float(row[4]))

def loadTerminalMismatchEEnergy(no_dangle=False, zero_energies=False):
    global terminalMismatchEEnergy
    terminalMismatchEEnergy = {}
    if no_dangle:
        filename = "inf_tstacke.csv"
    elif zero_energies:
        filename = "zero_tstacke.csv"
    else:
        filename = "tstacke.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                terminalMismatchEEnergy[str(row[0])+str(row[1])+str(row[2])+str(row[3])] = math.exp(-beta*float(row[4]))
                
def loadTerminalMismatchMEnergy(no_dangle=False, zero_energies=False):
    global terminalMismatchMEnergy
    terminalMismatchMEnergy = {}
    if no_dangle:
        filename = "inf_tstackm.csv"
    elif zero_energies:
        filename = "zero_tstackm.csv"
    else:
        filename = "tstackm.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                terminalMismatchMEnergy[str(row[0])+str(row[1])+str(row[2])+str(row[3])] = math.exp(-beta*float(row[4]))

# Fill in the hairpin energy dict (length : energy)
def loadHairpinEnergy(zero_energies=False):
    global hairpinEnergy
    hairpinEnergy = {}
    if zero_energies:
        filename = "zero_hairpin.csv"
    else:
        filename = "hairpin.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                hairpinEnergy[int(row[0])] = math.exp(-beta*float(row[1]))

# Fill in the special loop dict (sequence : energy)
def loadSpecialLoops(zero_energies=False):
    global specialLoops
    specialLoops = {}
    # Triloops
    # These are using the 2004 data (not currently implemented)
    '''with open("./data/triloop.csv") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                specialLoops[str(row[0])] = float(row[1])
    '''
    # Tetraloops
    if zero_energies:
        filename = "zero_tloop.csv"
    else:
        filename = "tloop.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                specialLoops[str(row[0])] = math.exp(-beta*float(row[1]))
    # Hexaloops
    # This actually does not belong here, but should be checked separately,
    # as hexaloops are from 2004 and have their final energy (not to be added to the earlier energy I don't *think*)
    '''with open("./data/hexaloop.csv") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the header row
            if i == 0:
                continue
            else:
                specialLoops[str(row[0])] = float(row[1])'''
# Fill in the C-loop dictionary (hairpinsequence : energy)
def loadCLoops(zero_energies=False):
    global CLoops
    CLoops = {}
    # C-only loop of size 3 (from Table 6 of Mathews)
    CLoops["CCC"] = math.exp(-beta*1.4)
    # C-only loops of larger size (up to 30) (from Table 6 of Mathews)
    if zero_energies:
        for i in CLoops:
            CLoops[i] = 1
    else:
        for i in range(4, 30):
            CLoops["C"*i] = math.exp(-beta*(0.3*i + 1.6))

# Asymmetry term
# Coef rounded from 0.48 to 0.5 to follow unafold and RNAeval with 99 rules
def loadAsymmetry(zero_energies=False):
    global asymmetry
    asymmetry = {}

    if zero_energies:
        for i in range(0, 31):
            asymmetry[i] = 1
    else:
        for i in range(0, 31):
            asymmetry[i] = math.exp(-beta*min(3, 0.5*abs(i)))
    #print("asymmetry", asymmetry)

# Fill the 1x1 internal loop dict (sequence : energy)
def loadInt11(zero_energies=False):
    global int11
    int11 = {}
    if zero_energies:
        filename = "zero_int11.csv"
    else:
        filename = "int11.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the first row
            if i == 0:
                continue
            else:
                int11[str(row[0])+str(row[1])+str(row[2])+
                      str(row[3])+str(row[4])+str(row[5])] = math.exp(-beta*float(row[6]))

# Fill the internal 2x1 bulge table (sequence : energy)
def loadInt21(zero_energies=False):
    global int21
    int21 = {}
    if zero_energies:
        filename = "zero_int21.csv"
    else:
        filename = "int21.csv"
    with open("./data/{}".format(filename)) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            # Skip the first row
            if i == 0:
                continue
            else:
                int21[str(row[0])+str(row[1])+str(row[2])+
                      str(row[3])+str(row[4])+str(row[5])+str(row[6])] = math.exp(-beta*float(row[7]))

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
                      str(row[4])+str(row[5])+str(row[6])+str(row[7])] = math.exp(-beta*float(row[8]))

# Fill the penalty dictionaries
def load_penalty_dicts(zero_energies=False):
    # Penalty for terminating a stack with anything but a GC
    # Hairpin version http://rna.urmc.rochester.edu/NNDB/turner99/hairpin-example-3.html
    # Rounded from 0.45 to 0.5 to follow unafold3.9
    # Switching it to 0.5 for testing against UNAFold
    global endPenalty
    # Loop version http://rna.urmc.rochester.edu/NNDB/turner99/internal-example-3.html
    # Rounded from 0.65 to 0.7 to follow unafold3.9 (tstacki.DAT) and RNAeval with 99 rules
    global loopEndPenalty
    # Loop version
    global loopFirstMismatchEnergy
    # Bonus for GU closure (closing pair with earlier two 5' bases before it : energy)
    # This is applied 'only to hairpins with a 5' closing G that is preceded by two G residues'
    # - Mathews
    global guClosureEnergy

    if zero_energies:
        endPenalty = {('AU'): 1, ('UA'): 1, ('GU'): 1, ('UG'): 1}
        loopEndPenalty = {('AU'): 1, ('UA'): 1, ('GU'): 1, ('UG'): 1}
        loopFirstMismatchEnergy = {('UU'): 1, ('GA'): 1, ('AG'): 1}
        guClosureEnergy = {('GGGU'): 1}
    else:
        endPenalty = {('AU'): math.exp(-beta*0.5), ('UA'): math.exp(-beta*0.5), ('GU'): math.exp(-beta*0.5), ('UG'): math.exp(-beta*0.5)}
        loopEndPenalty = {('AU'): math.exp(-beta*0.7), ('UA'): math.exp(-beta*0.7), ('GU'): math.exp(-beta*0.7), ('UG'): math.exp(-beta*0.7)}
        loopFirstMismatchEnergy = {('UU'): math.exp(-beta*(-0.7)), ('GA'): math.exp(-beta*(-1.1)), ('AG'): math.exp(-beta*(-1.1))}
        guClosureEnergy = {('GGGU'): math.exp(-beta*(-2.2))}


####################
# Global Variables #
####################

###########
# Testing #
###########

# Show the energy contributions
global debug_energy
debug_energy = False

# Show the structure
global debug_structure
debug_structure = False


############################
# Sequence and Environment #
############################

# The sequence itself
global s

# The length of the sequence
global N


# The scale constant (expected value for one base)
global scale
#scale = -0.34
scale = 0

# The maximum allowed internal loop or hairpin size
global l
l = 30


###################
#     Arrays      #
###################

# The four arrays (now five because of the two Z's)

# 1D, step through to calculate overall partition function ##BETTER EXPLANATION NEEDED##
global Z3
global Z5

# 2D, forces i and j to be paired
global Zb

# 2D, forces 'structure' somewhere between i and j
global Z1

# 2D, forces multibranch somewhere between i and j
global Z2


#####################
#  Bases and Pairs  #
#####################

# Convert base letter to corresponding int
global base_dict
base_dict = {'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'U': 3, 'u': 3, 'T': 3, 't': 3}

# List of lists of allowed pairs.
# R[i] contains all bases i is allowed to pair with to the right.
global R
R = []

# L[j] contains all bases j is allowed to pair with to the left.
global L
L = []

# Will contain the list of all pairs (for one given microstate calculated by the traceback)
global pair_dict

# Will contain the list of all dangles
global dangle_list

# Will contain the next e index for each d
global e_index


###################
#  Energy Tables  #
###################

# List of allowed pairs, stored as tuples
# Assumes only capital letters and no T's
global allowedPairs
allowedPairs = {('A','U'), ('U','A'), ('G','C'), ('C','G'), ('G','U'), ('U','G')}

# Stack energy
# Arranged like (5primeOuter5primeInner3primeInner3primeOuter : Energy)
global stackEnergy

# Internal Loops
global int11
global int21
global int22

# Hairpin
# Normal hairpin loops (just count how many in the loop)
global hairpinEnergy

# Special loops for which we know the specific energy
global specialLoops

# C-only loops
global CLoops

# Asymmetry value for bulge loops
global Asymmetry

# Internal loops
global internalLoopEnergy

# Bulge Loops
global bulgeLoopEnergy

# Dangle penalties
# Of the form (i, j, dangle)
global dangle5Energy

# Of the form (dangle, i, j)
global dangle3Energy

# Terminal Mismatch
# These are arranged like (5primeouter5primeinner3primeouter3primeinner : Energy)
global terminalMismatchHEnergy
global terminalMismatchIEnergy
global terminalMismatchEEnergy
global terminalMismatchMEnergy


####################
# Energy Penalties #
####################

# Multibranch Penalties
def load_multibranch_penalties(zero_energies=False):
    global ea
    global eb
    global ec

    if zero_energies:
        ea = 1
        eb = 1
        ec = 1
    else:
        ea = math.exp(-beta*3.4)
        eb = math.exp(-beta*0.0)
        ec = math.exp(-beta*0.4)


################################
# Loading and Useful functions #
################################

# Load in a sequence
def loadSequence(seq):
    global s
    global olds
    global N
    global Z3
    global Z5
    global Zb
    global Z1
    global Z2
    global e_index
    olds = str(seq).upper()
    N = len(olds)
    s = olds*2 # we double the string so that we can correctly handle the wraparound with R and L
    Z3 = np.full(N, 1, dtype=np.float64)
    Z5 = np.full(N, 1, dtype=np.float64)
    Zb = np.full((N,N), 0, dtype=np.float64)
    Z1 = np.full((N,N), 0, dtype=np.float64)
    Z2 = np.full((N,N), 0, dtype=np.float64)
    
# Add indices of other bases which can pair with this one
# R contains potential pairs to the right of the base
# L contains potintial pairs to the left of the base
# If given a Python list of lists of pairing partners,
# it will use those instead of calculating all legal pairs
# Usage: if base 0 is allowed to pair with 5 and 6, and base
# 1 is allowed to pair with base 7, it will look like:
# [[5,6], [7], ... ]
def fillAllowedPairsList(likely_pairs=None):
    global R
    global L
    R = []
    L = []

    # If we gave a list of pairs, load it
    if not not likely_pairs:
        for i in range(N):
            R.append([])
            L.append([])
            for base in likely_pairs[i]:
                if base >= min(i+4, N) and base < max(N, i+N-3):
                    R[-1].append(base)
                elif base+N >= min(i+4, N) and base < max(N, i+N-3): # for the wraparound case
                    R[-1].append(base+N)
            for base in likely_pairs[i]:
                if base >= 0 and base < i-3:
                    L[-1].append(base)
            L[-1].reverse()

    # Otherwise, load all legal pairs
    else:
        for i in range(N):
            # Make an empty list on the end of our list of lists
            R.append([])
            L.append([])
            
            # Fill it with the allowed pairs
            # These need to be AU, GC, or GU (or mirror/T-equivalent of those),
            # and must be able to form a helix of size at least 3
            for j in range(min(i+4, N), max(N, i+N-3)):
                if ( (s[i] in base_dict and s[j] in base_dict) and
                     (base_dict[s[i]] + base_dict[s[j]] == 3 or base_dict[s[i]] + base_dict[s[j]] == 5) ):
                    R[-1].append(j)

            for j in range(0, i-3):
                if ( (s[i] in base_dict and s[j] in base_dict) and 
                     (base_dict[s[i]] + base_dict[s[j]] == 3 or base_dict[s[i]] + base_dict[s[j]] == 5) ):
                    L[-1].append(j)
            L[-1].reverse()

# Fill the e_index matrix which gives the first eIdx given d and j
def filleIndex():
    global e_index
    e_index = {}
    #for d in range(1, N-5): # d can't be 0 because that's what i is.
                             # Similarly, it cannot go too far to the
                             # right or there will be no space for e and j.
    for d in range(-1, N): # it can go all the way to the end because of the wraparonud case
        if not R[d]: continue
        #print("d is ", d)
        eIdx = 0
        first = True
        #for j in range(d+5, N):
        for j in range(d+2, 2*N):
            #print("j is ", j)
            if first:
                if R[d][eIdx] < j:
                    first = False
                    e_index[(d, j)] = eIdx
            else:
                # make sure we're not off the end and step along by one if useful
                if eIdx + 1 < len(R[d]) and R[d][eIdx + 1] < j:
                    eIdx += 1
                e_index[(d, j)] = eIdx
        first = True

##################
#  Z1 Functions  #
##################

################
# Forward Fill #
################

# Note that i and j are NOT necessarily paired with one another
# We know there must be at least one pairing at some point between i and j though.

# We will make this faster by combining all the 'for' loops into one.
# So far have attempted to maintain readability by not doing this

# Get the 'only one stem' term
# Note: we multiply everything by ec at the end, which is why we don't apply it to each term
def singleStemTermF(i, j):
    # energy_sum will hold the sum over all k's. Then we'll multiply it by the Boltzmann factor.
    energy_sum = 0
    
    # No dangle
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k > min(j, N-1): break
        scale_num = j-k
        # Add in the paired value times the Boltzmann weight for the unpaired rest
        energy_sum += (Zb[i, k]
                      * (eb**(j-k))
                      * endPen)
    
    # 3'D dangle
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k > min(j-1, N-2): break
        scale_num = j-k
        energy_sum += (Zb[i, k]
                      * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                      * (eb**(j-k))
                      * endPen)
    # 5'D dangle
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k > min(j, N-1): break
        scale_num = j-k+1
        energy_sum += (dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                      * Zb[i+1, k]
                      * (eb**(j-k+1))
                      * endPen)
    # terminal dangle
    ## WE SHOULD QUESTION HERE WHETHER WE NEED TO CHECK IF THE DD COULD FORM A LEGAL PAIR, IN WHICH CASE WE DON'T ALLOW IT ##
    # or whether we check if it's in WC pairs #
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k > min(j-1, N-2): break
        scale_num = j-k+1
        # Here we treat the double dangle energy as a multibranch terminal mismatch energy ('tstackm')
        # ^^ We will check if this is correct
        if not (s[i], s[k+1]) in allowedPairs:
            energy_sum += (terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                          * Zb[i+1, k]
                          * (eb**(j-k+1))
                          * endPen)
        
    # Return the final Boltzmann factor (single scale factor for the single bonus)
    return ec*energy_sum

## DO WE NEED TO WORRY ABOUT THE CASE WHERE i==j ?? ##

# Get the 'more than one stem' term
def additionalStemTermF(i, j):
    # energy_sum will hold the sum of all k's. Then we'll multiply it by the Boltzmann factor, as before.
    energy_sum = 0
    
    # We loop over all internal k's, ranging from k = i to to k = j-1
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k > min(j-1, N-2): break
        # Add in the paired value times the 'some structure' value.
        energy_sum += (Zb[i, k]
                      * Z1[k+1, j]
                      * endPen)
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k > min(j-2, N-3): break
        energy_sum += (eb
                      * Zb[i, k]
                      * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                      * Z1[k+2, j]
                      * endPen)
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k > min(j-1, N-2): break
        energy_sum += (eb
                      * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                      * Zb[i+1, k]
                      * Z1[k+1, j]
                      * endPen)
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k > min(j-2, N-3): break
        if not (s[i], s[k+1]) in allowedPairs:
            energy_sum += (eb*eb
                          * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                          * Zb[i+1, k]
                          * Z1[k+2, j]
                          * endPen)

    # Return the final Biltzmann factor (single scale factor for the single bonus)
    return ec*energy_sum

# Get the 'i unpaired' term
def unpairedTerm1F(i, j):
    return eb*Z1[i+1, j]


##############
# Wraparound #
##############

def singleStemTermW(i, j):

    if i > N:
        raise Exception("something went wrong, i too big")

    # Base case
    if i == N:
        return 0

    # energy_sum will hold the sum over all k's. Then we'll multiply it by the Boltzmann factor.
    energy_sum = 0
   
    # No dangle
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k < N: continue
        if k > j: break
        scale_num = j-k
        # Add in the paired value times the Boltzmann weight for the unpaired rest
        energy_sum += (Zb[i, k-N]
                      * (eb**(j-k))
                      * endPen)

    # 3'D dangle
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k < N: continue
        if k > j-1: break
        scale_num = j-k
        energy_sum += (Zb[i, k-N]
                      * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                      * (eb**(j-k))
                      * endPen)
    # Base Case
    if i == N-1:
        return ec*energy_sum

    # 5'D dangle
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k < N: continue
        if k > j: break
        scale_num = j-k+1
        energy_sum += (dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                      * Zb[i+1, k-N]
                      * (eb**(j-k+1))
                      * endPen)
    # terminal dangle
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k < N: continue
        if k > j-1: break
        scale_num = j-k+1
        if not (s[i], s[k+1]) in allowedPairs:
            energy_sum += (terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                          * Zb[i+1, k-N]
                          * (eb**(j-k+1))
                          * endPen)
        
    # Return the final Boltzmann factor (single scale factor for the single bonus)
    return ec*energy_sum

# Get the 'more than one stem' term
def additionalStemTermW(i, j):
    #if i == 7 and j == 14:
        #print("got here")

    if i > N:
        raise Exception("something went wrong, i too big")

    # Base case
    if i == N:
        return 0

    # energy_sum will hold the sum of all k's. Then we'll multiply it by the Boltzmann factor, as before.
    energy_sum = 0
    
    # We loop over all internal k's, ranging from k = i to to k = j-1
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k == N-1: continue
        #if k > j-5: break replaced on 19 Jan
        if k >= N and k > j-5: break
        #if i == 7 and j == 14:
            #print("got here 2")
        # Add in the paired value times the 'some structure' value.
        energy_sum += (Zb[i, k-N]
                      * Z1[(k+1)%N, j-N]
                      * endPen)
        
    for k in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        if k == N-1 or k == N-2: continue
        #if k > j-6: break replaced on 19 Jan
        if k >= N and k > j-6: break
        #if i == 7 and j == 14:
            #print("got here 3")
        energy_sum += (eb
                      * Zb[i, k-N]
                      * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                      * Z1[(k+2)%N, j-N]
                      * endPen)

    # Base Case
    if i == N-1:
        return ec*energy_sum

    #print("i is", i)
    #print("s is", s)
    #print("len is", len(s))
    #print("len R", len(R))
    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k == N-1: continue
        #if k > j-5: break
        if k >= N and k > j-5: break
        #print("i, k, N are", i, k, N)
        #if i == 6 and j == 14:
            #print("got here 4")
        energy_sum += (eb
                      * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                      * Zb[i+1, k-N]
                      * Z1[(k+1)%N, j-N]
                      * endPen)

    for k in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        if k == N-1 or k == N-2: continue
        #if k > j-6: break
        if k >= N and k > j-6: break
        if not (s[i], s[k+1]) in allowedPairs:
            #if i == 7 and j == 14:
                #print("got here 5")
            energy_sum += (eb*eb
                          * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                          * Zb[i+1, k-N]
                          * Z1[(k+2)%N, j-N]
                          * endPen)

    # Return the final Biltzmann factor (single scale factor for the single bonus)
    #if i == 7 and j == 14: print("Z2 returning", energy_sum)
    #if i == 6 and j == 14: print("Z2 returning1", energy_sum)
    return ec*energy_sum

# Get the 'i unpaired' term
def unpairedTerm1W(i, j):

    if i > N:
        raise Exception("something went wrong, i too big")

    # Base case
    if i == N:
        return 0

    return eb*Z1[i+1, j-N]

##################
#  Z2 Functions  #
##################

################
# Forward Fill #
################

# The case in which i is unpaired (the earlier cases are already calculated in Z1)
def unpairedTerm2F(i, j):
    return eb*Z2[i+1, j]

##############
# Wraparound #
##############

# The case in which i is unpaired (the earlier cases are already calculated in Z1)
def unpairedTerm2W(i, j):

    if i > N:
        raise Exception("something went wrong, i too big")

    # Base case
    if i == N:
        return 0

    return eb*Z2[i+1, j-N]

# More Functions

# Fill the four arrays
# We do this by stepping along the diagonals, starting from the main diagonal
# (we don't worry about illegal hairpins because that is accounted for)

def fillMatrices1():
    # These indices make sure we're on the correct diagonals
    # (we don't worry about illegal hairpins because that is accounted for)
    # d is the difference between j and i
    for d in range(4, N):
        for i in range(N-d):
            j = i+d

            ########
            #  Zb  #
            ########
            
            # We fill up the first diagonal, and then every diagonal after it
            # This checks to make sure that the bases can pair (correct letters)
            #print()
            #if i == 0 and j == N-1:
                #print("Looking at %s,%s: %s"%(i, j, '('+s+')'))
            #else:
                #print("Looking at %s,%s: %s"%(i, j, s[:i]+'('+s[i:j+1]+')'+s[j+1:]))
                
            if not j in R[i]:
                Zb[i, j] = 0
            else:
                hTerm = hairpinTerm(i, j)
                #if hTerm > 0: print("hTerm {},{}:".format(i,j), hTerm)
                sTerm = stackTermF(i, j)
                #if sTerm > 0: print("sTerm {},{}:".format(i,j), sTerm)
                iTerm = internalTermF(i, j)
                #if iTerm > 0: print("iTerm {},{}:".format(i,j), iTerm)
                mTerm = multibranchTermF(i, j)
                #if mTerm > 0: print("mTerm {},{}:".format(i,j), mTerm)

                Zb[i, j] = hTerm + sTerm + iTerm + mTerm
                #if Zb[i, j] < 0:
                    #print("Hmm, less than 0")
                    #print("hairpinterm:  %s"%hTerm)
                    #print("stackterm:    %s"%sTerm)
                    #print("internalterm: %s"%iTerm)
                    #print("multiTerm:    %s"%mTerm)
                    #print("Total:        %s"%Zb[i, j])
                    #print()
                
            ########
            #  Z1  #
            ########
            
            # Calculate the additional stem term to save calculating it again for Z2
            additionalStem = additionalStemTermF(i, j)
            
            # Fill the diagonal
            Z1[i, j] = singleStemTermF(i, j) + additionalStem + unpairedTerm1F(i, j)
            
            ########
            #  Z2  #
            ########
            
            Z2[i, j] = additionalStem + unpairedTerm2F(i, j)


def fillMatrices2():
    # BASE CASES

    # These indices make sure we're on the correct diagonals
    # d is the difference between j and i
    #count = 1
    for d in reversed(range(1, N)):
        for i in range(d, N):
            j = N+i-d
            #print("i, j", i, j)
            #Zb[i, j-N] = count
            #count += 1
            #print("@ {}, {}".format(i, j))

            ########
            #  Zb  #
            ########
            # BASE CASES ARE HANDLED WITHIN THE INDIVIDUAL COMPONENT FUNCTIONS #
            
            # We fill up the first diagonal, and then every diagonal after it
            # This checks to make sure that the bases can pair (correct letters)
            #print()
            #if i == 0 and j == N-1:
                #print("Looking at %s,%s: %s"%(i, j, '('+s+')'))
            #else:
                #print("Looking at %s,%s: %s"%(i, j, s[:i]+'('+s[i:j+1]+')'+s[j+1:]))
            
            #print("i is", i)
            #print(R[i])
            if (j not in R[i]) or j >= i+N-3: # because then the hairpin would be too small
                #print("{} and {} with sep {} not good".format(s[i], s[j], i-j+N))
                Zb[i, j-N] = 0
                #print("{}, {} can't pair".format(i, j-N))
                # this could also be replaced with a "do nothing"
            else:
                eTerm = externalLoopTerm(i, j)
                #if i == 11 and j == 21: print("eTerm {},{}:".format(i,j), eTerm)
                sTerm = stackTermW(i, j)
                #if i == 11 and j == 21: print("sTerm {},{}:".format(i,j), sTerm)
                iTerm = internalTermW(i, j)
                #if i == 11 and j == 21: print("iTerm {},{}:".format(i,j), iTerm)
                mTerm = multibranchTermW(i, j)
                #if i == 11 and j == 21: print("mTerm {},{}:".format(i,j), mTerm)

                Zb[i, j-N] = eTerm + sTerm + iTerm + mTerm
                #print("Zb({},{})".format(i, j-N), Zb[i, j-N])
                #if Zb[i, j] < 0:
                    #print("Hmm, less than 0")
                    #print("hairpinterm:  %s"%hTerm)
                    #print("stackterm:    %s"%sTerm)
                    #print("internalterm: %s"%iTerm)
                    #print("multiTerm:    %s"%mTerm)
                    #print("Total:        %s"%Zb[i, j])
                    #print()
                
            ########
            #  Z1  #
            ########
            
            # Calculate the additional stem term to save calculating it again for Z2
            additionalStem = additionalStemTermW(i, j)
            
            # Fill the diagonal
            # Base Case
            if i == N-1:
                Z1[i, j-N] = singleStemTermW(i, j) + additionalStem

            else:
                Z1[i, j-N] = singleStemTermW(i, j) + additionalStem + unpairedTerm1W(i, j)
            
            ########
            #  Z2  #
            ########
            
            # Base Case
            if i == N-1:
                Z2[i, j-N] = additionalStem
            else:
                Z2[i, j-N] = additionalStem + unpairedTerm2W(i, j)
            
#################
# Z3 Functions  #
#################

# Fill the 1D Z3 matrix. We start at N and work our way down to the start.
## WE WILL NEED TO INCLUDE THE SCALE FACTOR CALCULATIONS HERE TOO AT SOME POINT? ##
def fillZ3():
    #print()
    #print("Now we fill the Z array, from right to left")
    global Z3
    # for i=N-5, N-6, ..., 0
    for i in reversed(range(N-4)):

        #Add in the AU penalty in the cases where we have a Zb
        energy_sum = 0

        # No dangles
        for k in R[i]:
            if k > N-3: break
            #print("  Looking at {}: {}".format(k, s[:i]+'('+'['+s[i]+']'+s[i+1:k]+'['+s[k]+']'+s[k+1:]+')'))
            endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
            energy_sum += Zb[i, k]*Z3[k+1]*endPen
            #print("  multiplying:")
            #try:
                #print("  Zb[i, k]: {}".format(-1/beta*math.log(Zb[i, k])))
            #except:
                #print("not allowed to pair")
            #print("  Z[k+1]:   {}".format(-1/beta*math.log(Z3[k+1])))
            #print("  EndPen:   {}".format(-1/beta*math.log(endPen)))

        # 3'D
        for k in R[i]:
            if k > N-3: break
            #print("  Looking at {}: {}".format(k, s[:i]+'('+'['+s[i]+']'+s[i+1:k]+'['+s[k]+']'+s[k+1:]+')'))
            # endPen is still for pair (i,k), so we don't need to recalculate it
            endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
            energy_sum += Zb[i, k]*dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]*Z3[k+2]*endPen
            #print("  multiplying:")
            #try:
                #print("  Zb[i, k]: {}".format(-1/beta*math.log(Zb[i, k])))
            #except:
                #print("not allowed to pair")
            #print("  3'dangle: {}".format(dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]))
            #print("  Z[k+2]:   {}".format(-1/beta*math.log(Z3[k+2])))
            #print("  EndPen:   {}".format(-1/beta*math.log(endPen)))

        # 5'D
        for k in R[i+1]:
            if k > N-3: break
            #print("  Looking at {}: {}".format(k, s[:i]+'('+s[i]+'['+s[i+1]+']'+s[i+2:k]+'['+s[k]+']'+s[k+1:]+')'))
            # Now looking at (i+1,k), so we calculate the new endPenalty
            endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
            energy_sum += Zb[i+1, k]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]*Z3[k+1]*endPen
            #print("  multiplying:")
            #try:
                #print("  Zb[i+1, k]: {}".format(-1/beta*math.log(Zb[i+1, k])))
            #except:
                #print("not allowed to pair")
            #print("  5'dangle: {}".format(dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]))
            #print("  Z[k+1]:   {}".format(-1/beta*math.log(Z3[k+1])))
            #print("  EndPen:   {}".format(-1/beta*math.log(endPen)))

        # ts
        # Note that if the terminal stack occurs with bases that can pair, we don't consider it
        for k in R[i+1]:
            if k > N-3: break
            #print("  Looking at {}: {}".format(k, s[:i]+'('+s[i]+'['+s[i+1]+']'+s[i+2:k]+'['+s[k]+']'+s[k+1:]+')'))
            #if str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1]) in terminalMismatchEEnergy: changed on 10 Jan
            # this is important because terminalMismatchEEnergy INCLUDES GU as a mismatch, which we don't allow
            if not (str(s[i]), str(s[k+1])) in allowedPairs:
                # ^^ we could replace this with simply checking if (s[k+1], s[i]) is a WC pair. We can check in the 'allowedPairs' list.
                # ^^ done on Jan 10, 2017
                endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
                energy_sum += (Zb[i+1, k]
                              * terminalMismatchEEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                              * Z3[k+2]
                              * endPen)
            else:
                energy_sum += 0 # If i+1 and k-1 can pair, then double dangle is not allowed
            #try:
                #print("   Zb[i+1, k]: {}".format(-1/beta*math.log(Zb[i+1, k])))
            #except:
                #print("not allowed to pair")
            #print("   Z[k+2]:     {}".format(-1/beta*math.log(Z3[k+2])))
            #print("   Ddangle:    {}".format(terminalMismatchEEnergy.get(str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1]), "Not Allowed")))
            #print("   EndPen:     {}".format(-1/beta*math.log(endPen)))

        # Case where k = N-2
        # No dangle and 3'D (we need to separate these because the endPen's change
        # shouldn't ^^^ this be a 5'D? (3 Jan)
        # ^ this is wrong, original is right (10 Jan)
        if N-2 in R[i]:
            # no dangle
            endPen = endPenalty.get(str(s[i])+str(s[N-2]), 1)
            energy_sum += Zb[i, N-2]*endPen
            #print("  Zb[i, N-2]:  {}".format(-1/beta*math.log(Zb[i, N-2])))
            #print("  EndPen:      {}".format(-1/beta*math.log(endPen)))

            # 3'D
            energy_sum += Zb[i, N-2]*dangle3Energy[str(s[i])+str(s[N-2])+str(s[N-1])]*endPen
            #print("  Zb[i, N-2]:  {}".format(-1/beta*math.log(Zb[i, N-2])))
            #print("  3dangle:     {}".format(-1/beta*math.log(dangle3Energy[str(s[i])+str(s[N-2])+str(s[N-1])])))
            #print("  EndPen:      {}".format(-1/beta*math.log(endPen)))

        if N-2 in R[i+1]:
            # 5'D
            endPen = endPenalty.get(str(s[i+1])+str(s[N-2]), 1)
            energy_sum += Zb[i+1, N-2]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-2])]*endPen
            #print("  Zb[i+1, N-2]:  {}".format(-1/beta*math.log(Zb[i+1, N-2])))
            #print("  5dangle:       {}".format(-1/beta*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-2])])))
            #print("  EndPen:        {}".format(-1/beta*math.log(endPen)))

            # tstacke
            #if str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1]) in terminalMismatchEEnergy: # changed on 10 Jan
            # important because terminalMismatchEEnergy contains GU pairs, which we don't want to consider
            if not (str(s[i]), str(s[N-1])) in allowedPairs:
                # ^^ we could replace this with simply checking if (s[k+1], s[i]) is a WC pair. We can check in the 'allowedPairs' list.
                # ^^ done
                # Use the same endPenalty as before, because we are still looking at (i, N-2)
                energy_sum += Zb[i+1, N-2]*terminalMismatchEEnergy[str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1])]*endPen
                #print("  Zb[i+1, N-2]:  {}".format(-1/beta*math.log(Zb[i+1, N-2])))
                #print("  tstack:        {}".format(-1/beta*math.log(terminalMismatchEEnergy[str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1])])))
                #print("  EndPen:        {}".format(-1/beta*math.log(endPen)))

        # Case where k = N-1 -- here we only have the no dangle and the 5'D cases
        if N-1 in R[i]:
            # No dangle
            endPen = endPenalty.get(str(s[i])+str(s[N-1]), 1)
            energy_sum += Zb[i, N-1]*endPen
            #print("  Zb[i, N-1]:  {}".format(Zb[i, N-1]))
            #print("  EndPen:      {}".format(-1/beta*math.log(endPen)))

        if N-1 in R[i+1]:
            # 5'D
            endPen = endPenalty.get(str(s[i+1])+str(s[N-1]), 1) # changed N-2 to N-1 on Jan 10. this was a bug.
            energy_sum += Zb[i+1, N-1]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-1])]*endPen
            #print("  Zb[i+1, N-1]:  {}".format(-1/beta*math.log(Zb[i+1, N-1])))
            #print("  5dangle:       {}".format(-1/beta*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-1])])))
            #print("  EndPen:        {}".format(-1/beta*math.log(endPen)))

        # Case where nothing is paired
        # Step along
        energy_sum += Z3[i+1]
    
        # Fill the array with the final value
        Z3[i] = energy_sum

#################
# Z5 Functions  #
#################

# Fill the 1D Z5 matrix. We start at 0 and work our way up to the end.
## WE WILL NEED TO INCLUDE THE SCALE FACTOR CALCULATIONS HERE TOO AT SOME POINT? ##
def fillZ5():
    #print()
    #print("Now we fill the Z5 array, from left to right")
    global Z5
    # for i = 4, 5, ... , N-1
    for j in range(4, N):
        #print(j)
        energy_sum = 0

        # No dangles
        #print(L[j])
        for k in L[j]:
            if k < 2: break
            endPen = endPenalty.get(str(s[k])+str(s[j]), 1)
            energy_sum += Zb[k, j]*Z5[k-1]*endPen
        # 3'D
        for k in L[j-1]:
            if k < 2: break
            endPen = endPenalty.get(str(s[k])+str(s[j-1]), 1)
            energy_sum += Zb[k, j-1]*dangle3Energy[str(s[k])+str(s[j-1])+str(s[j])]*Z5[k-1]*endPen
        # 5'D
        for k in L[j]:
            if k < 2: break
            endPen = endPenalty.get(str(s[k])+str(s[j]), 1)
            energy_sum += Zb[k, j]*dangle5Energy[str(s[k-1])+str(s[k])+str(s[j])]*Z5[k-2]*endPen
        # ts
        for k in L[j-1]:
            if k < 2: break
            #if str(s[j-1])+str(s[j])+str(s[k-1])+str(s[k]) in terminalMismatchEEnergy:
            if not (str(s[k-1]), str(s[j])) in allowedPairs: # done on Jan 10
                endPen = endPenalty.get(str(s[k])+str(s[j-1]), 1)
                energy_sum += (Zb[k, j-1]
                              * terminalMismatchEEnergy[str(s[j-1])+str(s[j])+str(s[k-1])+str(s[k])]
                              * Z5[k-2]
                              * endPen)
            else:
                energy_sum += 0

        # Case where k = 1
        # No dangle and 5'D (we need to separate these because the endPen's change
        if 1 in L[j]:
            endPen = endPenalty.get(str(s[1])+str(s[j]), 1)
            energy_sum += Zb[1, j]*endPen + Zb[1, j]*dangle5Energy[str(s[0])+str(s[1])+str(s[j])]*endPen

        # 3'D
        if 1 in L[j-1]:
            endPen = endPenalty.get(str(s[1])+str(s[j-1]), 1)
            energy_sum += Zb[1, j-1]*dangle3Energy[str(s[1])+str(s[j-1])+str(s[j])]*endPen 
            # tstacke
            #if str(s[j-1])+str(s[j])+str(s[0])+str(s[1]) in terminalMismatchEEnergy:
            if not (str(s[0]), str(s[j])) in allowedPairs: # done on Jan 10
                # Use the same endPenalty as before, because we are still looking at (1, j-1)
                energy_sum += Zb[1, j-1]*terminalMismatchEEnergy[str(s[j-1])+str(s[j])+str(s[0])+str(s[1])]*endPen

        # Case where k = 0 -- here we only have the no dangle and the 3'D cases
        # No dangle
        if 0 in L[j]:
            endPen = endPenalty.get(str(s[0])+str(s[j]), 1)
            energy_sum += Zb[0, j]*endPen

        # 3'D
        if 0 in L[j-1]:
            endPen = endPenalty.get(str(s[0])+str(s[j-1]), 1)
            energy_sum += Zb[0, j-1]*dangle3Energy[str(s[0])+str(s[j-1])+str(s[j])]*endPen

        # Case where nothing is paired
        # Step along
        energy_sum += Z5[j-1]
        #print(" Step along Z[i+1]:")
        #print("  {}".format(math.exp(-beta*(-scale))*Z[i+1]))
    
        # Fill the array with the final value
        Z5[j] = energy_sum


###############
#  Traceback  #
###############

# Given we've filled all the matrices, we calculate one of the possible structures

# Traceback for the Z
def traceZ3(i):
    #print()
    yes = False
    #print("Trace Z(%s)"%i)
    global Z3
    
    # This is where we'll record the indices of the dangling bases
    global dangle_list
    
    # Make sure we are not off the strand
    if i > N:
        raise ValueError("{} is larger than {}!".format(i, N))
        return "Something went wrong"

    # The base case: we're too close to the edge to have structure
    if i >=  N-4:
        if yes: print("We're at the end; returning 1")
        return 1
    
    # Choose a random number between 0 and 1
    rand = np.random.rand()
    #rand = 0.999999999999999
    if yes: print("picked %s"%rand)
    
    # We have two possibilities: i unpaired, or i paired with something in R[i]
    # Step
    prob0 = Z3[i+1]/Z3[i]
    if yes: print("Pr[Z3[{}]/Z3[{}]] = {}".format(i+1, i, prob0))
    
    # Check if it's big enough
    if yes: print(prob0, ">", rand, "?")
    if prob0 >= rand:
        if yes: print("Big enough!")
        if yes: print("Calling traceZ3({})".format(i+1))
        return traceZ3(i+1)
    if yes: print("Not enough, {} is paired".format(i))
    
    # If it's not, then look through the paired cases (ie the four dangle cases)
    # Same idea is fillZ3(), but we're looking at individual cases instead of adding all together.
    # (also everything is normalized by Z[i])
    #print(" No dangles:")
    if yes: print("Allowed Pairs: {}".format(R[i]))
    #print(R[1])
    for k in R[i]:
        if k > N-3: break
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        prob1 = Zb[i, k]*Z3[k+1]*endPen/Z3[i]
        if yes: print("Pr[(Zb[{}, {}]*Z[{}])/Z[{}]] = {}".format(i, k, k+1, i, prob1))
        if yes: print("Total: {}".format(prob0+prob1))
        
        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            if yes: print("Calling traceZb({}, {}) + traceZ({})".format(i, k, k+1))
            return traceZb(i, k) + traceZ3(k+1)
        #print("Not enough!")
        # Set prob0 to be the new baseline onto which we add
        prob0 = prob0 + prob1

    # 3'D
    #print(" 3'D (right):")
    #print("Allowed Pairs: {}".format(R[i]))
    for k in R[i]:
        if k > N-3: break # We're too close to the end
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        prob1 = Zb[i, k]*dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]*Z3[k+2]*endPen/Z3[i]
        if yes: print("Pr[(Zb[{}, {}]*Z[{}])/Z[{}]] = {}".format(i, k, k+2, i, prob1))
        if yes: print("Total: {}".format(prob0+prob1))
        
        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            # Add k+1 to the dangle list
            dangle_list.append(k)
            if yes: print("Added dangle {}".format(k+1))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {}) + traceZ({})".format(i, k, k+2))
            return traceZb(i, k) + traceZ3(k+2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # 5'D
    #print(" 5'D (left):")
    if yes: print("Allowed Pairs: {}".format(R[i+1]))
    for k in R[i+1]:
        if k > N-3: break
        endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
        prob1 = Zb[i+1, k]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]*Z3[k+1]*endPen/Z3[i]
        if yes: print("Pr[(Zb[{}, {}]*Z[{}])/Z[{}]] = {}".format(i+1, k, k+1, i, prob1))
        if yes: print("Total: {}".format(prob0+prob1))
        
        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            # Add i to the dangle list
            dangle_list.append(i)
            if yes: print("Added dangle {}".format(i))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {}) + traceZ({})".format(i+1, k, k+1))
            return traceZb(i+1, k) + traceZ3(k+1)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # ts
    # Note that if the double dangle occurs with bases that can pair, we don't consider it
    # (because then it counts as a pair, not a dangle)
    #print(" DD (both):")
    #print("Allowed Pairs: {}".format(R[i+1]))
    for k in R[i+1]:
        if k > N-3: break
        #if any(k+1 in pair for pair in pair_dict): break # We can't have a dangle if it's already paired
        #if str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1]) in terminalMismatchEEnergy:
        if not (str(s[i]), str(s[k+1])) in allowedPairs: # done on Jan 10
            endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1) # THIS WAS NOT HERE BEFORE, ADDED FEB 19
            prob1 = Zb[i+1, k]*terminalMismatchEEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]*Z3[k+2]*endPen/Z3[i]
        else:
            prob1 = 0 # If i+1 and k-1 can pair, then double dangle is not allowed
        if yes: print("Pr[(Zb[{}, {}]*Z[{}])/Z[{}]] = {}".format(i+1, k, k+2, i, prob1))
        if yes: print("Total: {}".format(prob0+prob1))
        
        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            # Add i and k+1 to the dangle list
            dangle_list.append(i)
            dangle_list.append(k)
            if yes: print("Added dangles {}, {}".format(i, k+1))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {}) + traceZ({})".format(i+1, k, k+2))
            return traceZb(i+1, k) + traceZ3(k+2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # Case where k = N-2
    # Note: from here on out, we do not return any traceZ3, because we are close enough to N-1
    # that it is not possible for it to yield another Zb.
    # No dangle and 3'D (we need to separate these because the endPen's change
    if N-2 in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[N-2]), 1)
        prob1 = Zb[i, N-2]*endPen/Z3[i]
         
        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            if yes: print("Calling traceZb({}, {})".format(i, N-2))
            return traceZb(i, N-2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

        prob1 = Zb[i, N-2]*dangle3Energy[str(s[i])+str(s[N-2])+str(s[N-1])]*endPen/Z3[i]

        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            dangle_list.append(N-2)
            if yes: print("Added dangle {}".format(N-2))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {})".format(i, N-2))
            return traceZb(i, N-2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # 5'D
    if N-2 in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[N-2]), 1)
        prob1 = Zb[i+1, N-2]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-2])]*endPen/Z3[i]

        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            dangle_list.append(i)
            if yes: print("Added dangle {}".format(i))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {})".format(i+1, N-2))
            return traceZb(i+1, N-2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

        # tstacke
        #if str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1]) in terminalMismatchEEnergy:
        if not (str(s[i]), str(s[N-1])) in allowedPairs: # done on Jan 10
            # ^^ we could replace this with simply checking if (s[k+1], s[i]) is a WC pair. We can check in the 'allowedPairs' list.
            # Use the same endPenalty as before, because we are still looking at (i, N-2)
            prob1 = Zb[i+1, N-2]*terminalMismatchEEnergy[str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1])]*endPen/Z3[i]
        else:
            prob1 = 0

        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            dangle_list.append(i)
            dangle_list.append(N-2)
            if yes: print("Added dangles {}, {}".format(i, N-2))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {})".format(i+1, N-2))
            return traceZb(i+1, N-2)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # Case where k = N-1 -- here we only have the no dangle and the 5'D cases
    # No dangle
    if N-1 in R[i]:
        endPen = endPenalty.get(str(s[i])+str(s[N-1]), 1)
        prob1 = Zb[i, N-1]*endPen/Z3[i]

        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            if yes: print("Calling traceZb({}, {})".format(i, N-1))
            return traceZb(i, N-1)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # 5'D
    if N-1 in R[i+1]:
        endPen = endPenalty.get(str(s[i+1])+str(s[N-1]), 1) # fixed on Jan 10 2017. this was a bug.
        prob1 = Zb[i+1, N-1]*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-1])]*endPen/Z3[i]

        # Check if it's big enough
        if prob0 + prob1 >= rand:
            if yes: print("Big enough!")
            dangle_list.append(i)
            if yes: print("Added dangle {}".format(i))
            if yes: print("dangle_list:", dangle_list)
            if yes: print("Calling traceZb({}, {})".format(i+1, N-1))
            return traceZb(i+1, N-1)
        if yes: print("Not enough!")
        prob0 = prob0 + prob1

    # If we get here, then something is wrong with our probabilities, 
    # which should have added to one and thus cover all potential values of rand
    #print("random number was %s"%rand)
    #print("we got up to      %s"%(prob0))
    raise Exception("Did not reach random number in probability sum")
    return "Something went wrong"

# Traceback for Zb
def traceZb(i, j):
    yes = False
    #print("Trace Zb({}, {})".format(i, j))
    # i,j is a pair. Add them to the pair list!
    global pair_dict
    global dangle_list
    if yes: print("dangle list", dangle_list)
    pair_dict[i] = j
    
    # Choose a random number between 0 and 1
    rand = np.random.rand()
    #rand = 0.999999999999999
    #print("rand = {}".format(rand))
    
    # We have a few possibilities. In order of expected likelyhood:
    # i, j stack;
    # i, j close hairpin;
    # i, j close internal loop;
    # i, j close multibranch loop.

    # Check the stacked case
    # Note that if it cannot stack this will be taken care of by the multiplication by Zb[i+1, j-1]
    # Stacking case
    prob0 = stackTermF(i, j)/Zb[i, j]
    #print("stackTermF({},{}/Zb[{},{}]) = {}".format(i, j, i, j, prob0))
    if prob0 >= rand:
        #print("Big enough!")
        #print("Calling traceZb({}, {})".format(i+1, j-1))
        return traceZb(i+1, j-1)
    #print("Not enough!")
        
    # Hairpin case
    prob1 = hairpinTerm(i, j)/Zb[i, j]
    #print("hairpinTerm({},{})/Zb[{},{}] = {}".format(i, j, i, j, prob1))
    #print("Total:   {}".format(prob0+prob1))
    if prob0+prob1 >= rand:
        #print("Big enough!")
        #print("Hairpin; calling nothing")
        # Do nothing
        return 0
    #print("Not enough!")
    prob0 = prob0 + prob1
    
    for d in range(i+1, min(i+l,j)):

        # We use this variable for the duration of this loop
        if not (d, j) in e_index: continue
        #print("d is ", d)
        #print("R[d] is ", R[d])
        eIdx = e_index[(d, j)]
        #print("starting eIdx is ", eIdx)
        e = R[d][eIdx]
        #print("starting e is ", e)
        
        while j-e+d-i < l+2 and eIdx >= 0:
            #print("eIdx is ", eIdx)
            #print("d, e are {}, {}".format(d, e))
            #print("we are at d = {}".format(d))
            #print("R[{}] = {}".format(d, R[d]))
            
            # We get the actual value of e
            e = R[d][eIdx]

            # Set the energy of this microstate (Jan 16)
            energyZ = 1 # as of feb 19 there is no energyG
            #print("new e is ", e)

            #print("e is %d"%e)
            if e >= j: break
                
            # Break if loop is too big
            if d-i-1 + j-e-1 >= l: break

            # The scale number is (d-i) + (j-e)
            scale_num = (d-i) + (j-e)
            # Number in top loop
            top = d-i-1
            # Number in bottom loop
            bottom = j-e-1
                        
            # Base situation (we don't want to count this)
            if top == 0 and bottom == 0:
                eIdx -= 1
                continue
            # Bulge situation
            elif top == 0 or bottom == 0:
                energyZ *= bulgeLoopEnergy[top+bottom]
                # End Penalty
                if top > 1 or bottom > 1:
                    energyZ *= endPenalty.get(str(s[i])+str(s[j]), 1)
                    energyZ *= endPenalty.get(str(s[e])+str(s[d]), 1)
                # Stacking energy (added 5 Dec 2016)
                if top == 1 or bottom == 1:
                    energyZ *= stackEnergy.get(str(s[i])+str(s[d])+str(s[e])+str(s[j]), 1)
            # 1x1 situation
            elif top == 1 and bottom == 1:
                energyZ *= int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])]
            # 2x1 situation
            elif top == 2 and bottom == 1:
                energyZ *= int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])]
            # 1x2 situation
            elif top == 1 and bottom == 2:
                energyZ *= int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
            # 2x2 situation
            elif top == 2 and bottom == 2:
                energyZ *= int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
            # General case
            # See http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
            else:
                # Initiation term
                energyZ *= internalLoopEnergy[top+bottom]
            
                # Asymmetry term
                # Coef rounded from 0.48 to 0.5 to follow unafold and RNAeval with 99 rules
                energyZ *= asymmetry[abs(top-bottom)]
            
                # 'First Mismatch' term
                if top > 1 and bottom > 1:
                    # 5' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 1)
                    # 3' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 1)
      
                # Loop closure penalty term
                energyZ *= loopEndPenalty.get(str(s[i])+str(s[j]), 1)
                energyZ *= loopEndPenalty.get(str(s[e])+str(s[d]), 1)
              
            # Turn the energy into Boltzmann form
            energyZ *= Zb[d, e]
            
            #print("energyZ = {}".format(energy_Z))

            # Now we check the probabilities
            prob1 = energyZ/Zb[i, j]
            #print("({},{},{},{}) loop = {}".format(i, d, e, j, prob1))
            #print("Total = {}".format(prob0+prob1))
            if prob0+prob1 >= rand:
                #print("Big enough!")
                #print("Calling traceZb({},{})".format(d,e))
                return traceZb(d, e)
            #print("Not enough")
            prob0 = prob0 + prob1
            
            # Change the e_index
            eIdx -= 1
            #print("new eIdx is", eIdx)
        
    # Multibranch case
    # Will add scale factor later

    # No dangle
    prob1 = ea*ec*Z2[i+1, j-1]/Zb[i, j]
    #print("Multibranch no dangle = {}".format(prob1))
    #print("Total: {}".format(prob0+prob1))
    if prob0+prob1 >= rand:
        #print("Big enough!")
        #print("Adding nothing to dangle list")
        #print("Calling traceZ2({},{})".format(i+1, j-1))
        return traceZ2(i+1, j-1)
    #print("Not enough")
    prob0 = prob0 + prob1
    
    # 3'D
    prob1 = ea*eb*ec*dangle3Energy[str(s[j])+str(s[i])+str(s[i+1])]*Z2[i+2, j-1]/Zb[i, j]
    #print("Multibranch with 3'D = {}".format(prob1))
    #print("Total: {}".format(prob0+prob1))
    if prob0+prob1 >= rand:
        #print("Big enough!")
        #print("Adding {} to dangle list".format(i+1))
        dangle_list.append(i)
        #print("Calling traceZ2({},{})".format(i+2, j-1))
        return traceZ2(i+2, j-1)
    #print("Not enough")
    prob0 = prob0 + prob1
    
    # 5'D
    prob1 = ea*eb*ec*dangle5Energy[str(s[j-1])+str(s[j])+str(s[i])]*Z2[i+1, j-2]/Zb[i, j]
    #print("Multibranch with 5'D = {}".format(prob1))
    #print("Total: {}".format(prob0+prob1))
    if prob0+prob1 >= rand:
        #print("Big enough!")
        #print("Adding {} to dangle list".format(j-1))
        dangle_list.append(j-1)
        #print("Calling traceZ2({},{})".format(i+1, j-2))
        return traceZ2(i+1, j-2)
    #print("Not enough")
    prob0 = prob0 + prob1
    
    # DD
    ## IF WE ARE CHECKING FOR ALLOWED PAIRS, WE MUST DO IT HERE TOO ##
    #if 1:#(s[i+1], s[j-1]) not in allowedPairs:
    if not (s[i+1], s[j-1]) in allowedPairs:
        prob1 = (ea*eb*eb*ec*terminalMismatchMEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]
                * Z2[i+2, j-2])/Zb[i, j]
        #print("Multibranch with DD = {}".format(prob1))
    #print("Total: {}".format(prob0+prob1))
    if prob0+prob1 >= rand:
        #print("Big enough!")
        #print("Adding {} and {} to the dangle list".format(i+1, j-1))
        dangle_list.append(i)
        dangle_list.append(j-1)
        #print("Calling traceZ2({},{})".format(i+2, j-2))
        return traceZ2(i+2, j-2)
    #print("Not enough")
    prob0 = prob0 + prob1
    
    # If we get here, then something is wrong with our probabilities, 
    # which should have added to one and thus cover all potential values of rand
    #print("random number was {}".format(rand))
    #print("we got up to      {}".format(prob0))
    raise Exception("Did not reach random number in probability sum")
    return "Something went wrong"

# Traceback for Z1
def traceZ1(i, j):

    global dangle_list

    # Base case
    if i == j: 
        #print("Base case, returning 0")
        return 0
    
    # Error case
    if i > j:
        raise ValueError("i greater than j, something is wrong")
        return "Something went wrong"
    
    # Choose a random number between 0 and 1
    rand = np.random.rand()
    #rand = 0.99999
    #print("We chose {}".format(rand))
    
    # We have a few possibilities. In order of expected likelyhood:
    # i unpaired    
    # i, k paired, later structure
    # i, k paired, no later structure

    # i unpaired
    prob0 = eb*Z1[i+1, j]/Z1[i, j]
    #print("Unpaired = {}".format(prob0))
    #print("Total = {}".format(prob0))
    if prob0 >= rand:
        #print("Big enough!")
        #print("Calling traceZ({})".format(i+1))
        return traceZ1(i+1, j)
    #print("Not enough")
    
    # i, k paired, no later structure
    #print("  No later structure")
    #print(" No dangle")
    for k in R[i]:
        #print("Allowed Pairs: {}".format(R[i]))
        if k > min(j, N-1): break
        # How many are we adding on?
        scale_num = j-k
        # Add in the paired value times the Boltzmann weight for the unpaired rest
        prob1 = ec*Zb[i, k]*(eb**(j-k))/Z1[i, j]
        #print("Zb[{},{}]*exp/Z1[{},{}] = {}".format(i, k, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Calling Zb({},{})".format(i, k))
            return traceZb(i, k)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print(" 3'D")
    for k in R[i]:
        if k > min(j-1, N-2): break
        scale_num = j-k
        prob1 = (ec
                * Zb[i, k]
                * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                * (eb**(j-k)))/Z1[i, j]
        #print("Zb[{},{}]*exp*3'D/Z1[{},{}] = {}".format(i, k, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(k+1))
            #print("Calling Zb({},{})".format(i, k))
            dangle_list.append(k)
            return traceZb(i, k)
        #print("Not enough")
        prob0 = prob0 + prob1
        
    #print(" 5'D")
    for k in R[i+1]:
        if k > min(j, N-1): break
        scale_num = j-k+1
        prob1 = (ec
                 *dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                 * Zb[i+1, k]
                 * (eb**(j-k+1)))/Z1[i, j]
        #print("Zb[{},{}]*exp*5'D/Z1[{},{}] = {}".format(i+1, k, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(i))
            #print("Calling Zb({},{})".format(i+1, k))
            dangle_list.append(i)
            return traceZb(i+1, k)
        #print("Not enough")
        prob0 = prob0 + prob1
        
    #print(" ts")
    for k in R[i+1]:
        if k > min(j-1, N-2): break
        scale_num = j-k+1
        # Here we treat the double dangle energy as a multibranch terminal mismatch energy ('tstackm')
        if not (s[i], s[k+1]) in allowedPairs:
            prob1 = (ec
                    * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                    * Zb[i+1, k]
                    * (eb**(j-k+1)))/Z1[i, j]
        else:
            prob1 = 0
            #print("Zb[{},{}]*exp*DD/Z1[{},{}] = {}".format(i+1, k, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} and {} to the dangle list".format(i, k+1))
            #print("Calling Zb({},{})".format(i+1, k))
            dangle_list.append(i)
            dangle_list.append(k)
            return traceZb(i+1, k)
        #print("Not enough")
        prob0 = prob0 + prob1

    # i, k paired, some later structure
    #print("  Some later structure")
    #print(" No dangle")
    for k in R[i]:
        if k > min(j-1, N-2): break
        # Add in the paired value times the 'some structure' value.
        prob1 = ec*Zb[i, k]*Z1[k+1, j]/Z1[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z1[{},{}] = {}".format(i, k, k+1, j, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Calling Zb({},{}) and Z1({},{})".format(i, k, k+1, j))
            return traceZb(i, k) + traceZ1(k+1, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print("3'D")
    for k in R[i]:
        if k > min(j-2, N-3): break
        prob1 = (eb*ec*Zb[i, k]
                * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                * Z1[k+2, j])/Z1[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z1[{},{}] = {}".format(i,k,k+2,j,i,j))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(k+1))
            #print("Calling Zb({},{}) and Z1({},{})".format(i, k, k+2, j))
            dangle_list.append(k)
            return traceZb(i, k) + traceZ1(k+2, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print("5'D")
    for k in R[i+1]:
        if k > min(j-1, N-2): break
        prob1 = (eb*ec
                * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                * Zb[i+1, k]
                * Z1[k+1, j])/Z1[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z1[{},{}]".format(i+1,k,k+1,j,i,j))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(i))
            #print("Calling Zb({},{}) and Z1({},{})".format(i+1, k, k+1, j))
            dangle_list.append(i)
            return traceZb(i+1, k) + traceZ1(k+1, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print("ts")
    for k in R[i+1]:
        if k > min(j-2, N-3): break
        if not (s[i], s[k+1]) in allowedPairs:
            prob1 = (eb*eb*ec
                    * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                    * Zb[i+1, k]
                    * Z1[k+2, j])/Z1[i, j]
        else:
            prob1 = 0
        #print("Zb[{},{}]*Z1[{},{}]/Z1[{},{}]".format(i+1,k,k+2,j,i,j))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} and {} to the dangle list".format(i, k+1))
            #print("Calling Zb({},{}) and Z1({},{})".format(i+1, k, k+2, j))
            dangle_list.append(i)
            dangle_list.append(k)
            return traceZb(i+1, k) + traceZ1(k+2, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    # If we get here, then something is wrong with our probabilities, 
    # which should have added to one and thus cover all potential values of rand
    #print("random number was %s"%rand)
    #print("we got up to      %s"%(prob0))
    raise Exception("Did not reach random number in probability sum")
    return "Something went wrong"

# Traceback for Z2
def traceZ2(i, j):

    global dangle_list

    # Base case
    if i == j: return 0
    
    # Error case
    if i > j:
        raise ValueError("i is greater than j, something is wrong")
        return "Something went wrong"
    
    # Choose a random number between 0 and 1
    rand = np.random.rand()
    #rand = 0.9999999
    #print("We chose {}".format(rand))
    
    # We have a two possibilities. In order of expected likelyhood:
    # i unpaired    
    # i, k paired, some later structure
    
    # i unpaired
    prob0 = eb*Z2[i+1, j]/Z2[i, j]
    #print("Unpaired: {}".format(prob0))
    #print("Total = {}".format(prob0))
    if prob0 >= rand:
        #print("Big enough!")
        #print("Calling traceZ2({}, {})".format(i+1, j))
        return traceZ2(i+1, j)
    
    # i, k paired, with later structure
    #print("  Some later structure")
    #print(" No dangle")
    for k in R[i]:
        if k > min(j-1, N-2): break
        # Add in the paired value times the 'some structure' value.
        prob1 = ec*Zb[i, k]*Z1[k+1, j]/Z2[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z2[{},{}] = {}".format(i, k, k+1, j, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Calling Zb({},{}) and Z1({},{})".format(i, k, k+1, j))
            return traceZb(i, k) + traceZ1(k+1, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print(" 3'D")
    for k in R[i]:
        if k > min(j-2, N-3): break
        prob1 = (eb*ec
                * Zb[i, k]
                * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                * Z1[k+2, j])/Z2[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z2[{},{}] = {}".format(i, k, k+2, j, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
                
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(k+1))
            #print("Calling Zb({},{}) and Z1({},{})".format(i, k, k+2, j))
            dangle_list.append(k)
            return traceZb(i, k) + traceZ1(k+2, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print(" 5'D")
    for k in R[i+1]:
        if k > min(j-1, N-2): break
        prob1 = (eb*ec
                * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                * Zb[i+1, k]
                * Z1[k+1, j])/Z2[i, j]
        #print("Zb[{},{}]*Z1[{},{}]/Z2[{},{}] = {}".format(i+1, k, k+1, j, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} to the dangle list".format(i))
            #print("Calling Zb({},{}) and Z1({},{})".format(i+1, k, k+1, j))
            dangle_list.append(i)
            return traceZb(i+1, k) + traceZ1(k+1, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    #print(" ts")
    for k in R[i+1]:
        if k > min(j-2, N-3): break
        if not (s[i], s[k+1]) in allowedPairs:
            prob1 = (eb*eb*ec
                    * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                    * Zb[i+1, k]
                    * Z1[k+2, j])/Z2[i, j]
        else:
            prob1 = 0
            #print("Zb[{},{}]*Z1[{},{}]/Z2[{},{}] = {}".format(i+1, k, k+2, j, i, j, prob1))
        #print("Total = {}".format(prob0+prob1))
        
        # Check if it's enough
        if prob0 + prob1 >= rand:
            #print("Big enough!")
            #print("Added {} and {} to the dangle list".format(i, k+1))
            #print("Calling Zb({},{}) and Z1({},{})".format(i+1, k, k+2, j))
            dangle_list.append(i)
            dangle_list.append(k)
            return traceZb(i+1, k) + traceZ1(k+2, j)
        #print("Not enough")
        prob0 = prob0 + prob1
    
    # If we get here, then something is wrong with our probabilities, 
    # which should have added to one and thus cover all potential values of rand
    #print("random number was %s"%rand)
    #print("we got up to      %s"%(prob0))
    raise Exception("Did not reach random number in probability sum")
    return "Something went wrong"

################
#  TracebackD  #
################

# The Deterministic Traceback. Given a list of pairs and a last of dangles,
# this traceback returns the partition function associated with that structure
# Traceback for the Z
def traceZ3d(i):
    if debug_structure: print("Z3d({})".format(i))

    # We'll read from the dangle and pair lists to see what needs to be paired/stacked
    # Then we multiply the pre-existing trace_energy by the resulting energy (in 
    # partition function form)
    global dangle_list
    global pair_dict
    
    # Make sure we are not off the strand
    if i > N:
        raise ValueError("{} is larger than {}!".format(i, N))
        return "Something went wrong"
    # The base case: we're too close to the edge to have structure
    if i >=  N-4:
        return 1
 
    # We have two possibilities: i unpaired, or i paired with something in R[i]

    # If this base is not in the list of paired bases
    # No change in energy due to this - just move along by 1
    if not i in pair_dict and not (i+1 in pair_dict and i in dangle_list):
        if debug_energy: print("Z3 steps past {}:".format(i), 1)
        return traceZ3d(i+1)
    
    # If i is paired, do this:
    if i in pair_dict:

        if i-1 in dangle_list:
            raise Exception("This should not have happened")
        k = pair_dict[i]

        # General Case
        # No dangle
        if k <= N-3:
            #print("general case")
            endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
            if not k in dangle_list:
                if debug_energy: print("Z3 endPen {},{}:".format(i, k), (-1/beta)*math.log(endPen))
                return traceZbd(i)*traceZ3d(k+1)*endPen

        # 3'D
            if k in dangle_list:
                if debug_energy: 
                    print("Z3 endPen {},{}:".format(i, k), (-1/beta)*math.log(endPen))
                    print("Z3 dangle3 {},{},{}:".format(i, k, k+1), (-1/beta)*math.log(dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]))
                return (traceZbd(i)
                       * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                       * traceZ3d(k+2)
                       * endPen)
        # NOTE: 5'D case handled later (when we set j = pair_dict[i+1])

        # k == N-2 Base Case
        # No dangle 
        if k == N-2:
            #print("k == N-2")
            endPen = endPenalty.get(str(s[i])+str(s[N-2]), 1)
            if not N-2 in dangle_list:
                if debug_energy:
                    print("Z3 endPen {},{}:".format(i, N-2), (-1/beta)*math.log(endPen))
                return traceZbd(i)*endPen

        # 3'D 
            if N-2 in dangle_list:
                if debug_energy:
                    print("Z3 endPen {},{}:".format(i, N-2), (-1/beta)*math.log(endPen))
                    print("Z3 dangle3 {},{},{}:".format(i, N-2, N-1), (-1/beta)*math.log(dangle3Energy[str(s[i])+str(s[N-2])+str(s[N-1])]))
                return traceZbd(i)*dangle3Energy[str(s[i])+str(s[N-2])+str(s[N-1])]*endPen
        
        # k == N-1 Base Case
        # No dangle
        if k == N-1:
            #print("k == N-1")
            endPen = endPenalty.get(str(s[i])+str(s[N-1]), 1)
            if not N-1 in dangle_list:
                if debug_energy: print("Z3 endPen {},{}:".format(i, N-1), (-1/beta)*math.log(endPen))
                return traceZbd(i)*endPen

    # Otherwise, i+1 must be paired with the new k
    else:
        k = pair_dict[i+1]
        if not i in dangle_list:
            raise Exception("It should have been here")

        # General Case
        if k <= N-3:
            endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
            # 5'D
            if not k in dangle_list:
                if debug_energy: 
                    print("Z3 endPen {},{}:".format(i+1, k), (-1/beta)*math.log(endPen))
                    print("Z3 dangle5 {},{},{}:".format(i, i+1, k), (-1/beta)*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]))
                return traceZbd(i+1)*dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]*traceZ3d(k+1)*endPen
            
            # ts
            if k in dangle_list:
                #if any(k+1 in pair for pair in pair_dict): break # We can't have a dangle if it's already paired
                #if str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1]) in allowedPairs: changed on 10 Jan
                # this was simply incorrect
                if (str(s[i]), str(s[k+1])) in allowedPairs:
                    raise Exception("that's not allowed")
                if debug_energy:
                    print("Z3 endPen {},{}:".format(i+1, k), (-1/beta)*math.log(endPen)) 
                    print("Z3 ts {},{},{},{}:".format(k, k+1, i, i+1), (-1/beta)*math.log(terminalMismatchEEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]))
                return traceZbd(i+1)*terminalMismatchEEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]*traceZ3d(k+2)*endPen
            
        # Base case k == N-2
        # 5'D
        if k == N-2:
            endPen = endPenalty.get(str(s[i+1])+str(s[N-2]), 1)
            if not N-2 in dangle_list:
                if debug_energy: 
                    print("Z3 endPen {},{}:".format(i+1, N-2), (-1/beta)*math.log(endPen))
                    print("Z3 dangle5 {},{},{}:".format(i, i+1, N-2), (-1/beta)*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-2])]))
                return traceZbd(i+1)*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-2])]*endPen

        # ts
            if N-2 in dangle_list:
                #if str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1]) in allowedPairs: # changed on 10 Jan because it didn't make sense
                if (str(s[i]), str(s[N-1])) in allowedPairs:
                    raise Exception("not allowed")
                if debug_energy:
                    print("Z3 endPen {},{}:".format(i+1, N-2), (-1/beta)*math.log(endPen))
                    print("Z3 ts {},{},{}:".format(N-2, N-1, i, i+1), (-1/beta)*math.log(terminalMismatchEEnergy[str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1])]))
                return traceZbd(i+1)*terminalMismatchEEnergy[str(s[N-2])+str(s[N-1])+str(s[i])+str(s[i+1])]*endPen

        # Base Case k == N-1
        # 5'D
        else:
            endPen = endPenalty.get(str(s[i+1])+str(s[N-1]), 1)
            if debug_energy: 
                print("Z3 endPen {},{}:".format(i, N-1), (-1/beta)*math.log(endPen))
                print("Z3 dangle5 {},{},{}:".format(i, i+1, N-1), (-1/beta)*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-1])]))
            #print("calling on", i+1, pair_dict[i+1])
            return traceZbd(i+1)*dangle5Energy[str(s[i])+str(s[i+1])+str(s[N-1])]*endPen

    # If we get here, then something is wrong with our dangle_list or pair_dict, 
    # which should have added to one and thus cover all potential values of rand
    print("sequence", olds)
    print("pair_dict", pair_dict)
    print("dangle_list", dangle_list)
    raise Exception("Did not reach random number in probability sum")
    return "Something went wrong"

# Given a place to start looking, returns the number of structures inside this one
# Can be arbitrarily large
def get_num_structures(i):
    closing_base = pair_dict[i]
    next_base = i+1
    num_structures = 0
    
    while next_base < closing_base:
        if next_base in pair_dict:
            num_structures += 1
            next_base = pair_dict[next_base] + 1
        else:
            next_base += 1
    #print("num_structures", num_structures)
    return num_structures

# Givena place to start looking, returns the index of the next starting loop
def get_next_loop(i):
    j = pair_dict[i]
    next_loop = i+1

    while next_loop < j-4:
        if next_loop in pair_dict:
            return next_loop
        else:
            next_loop += 1
    raise Exception("Didn't find next structure. This is unexpected")
    return -1

# Traceback for Zb
# Note that we don't need a j because we already know what it is 
# because j = pair_dict[i]
def traceZbd(i):
    global pair_dict
    global dangle_list

    if debug_structure: print("Zbd({}, {})".format(i, pair_dict[i]))

    j = pair_dict[i]
    if not (s[i], s[j]) in allowedPairs:
        print("pair_dict[{}]".format(i), j)
        print("sequence", olds)
        print("allowedPairs", allowedPairs)
        raise Exception

    # Which case are we in?
    # case = 0: i, j close hairpin;
    # case = 1: i, j close internal loop or stack;
    # case >= 2: i, j close multibranch loop.
    case = get_num_structures(i)

    # Hairpin loop case
    if case == 0:
        #print("case 0")
        #print("hairpin multiple of", hairpinTerm(i, j))
        if debug_energy: print("hairpin of:", (-1/beta) * math.log(hairpinTerm(i, j)))
        return hairpinTerm(i, j)
    
    # Stacking loop case
    if case == 1 and (i+1, j-1) in pair_dict.items():
        # We're adding on 2 bases
        scale_num = 2

        # Get the stacking energy
        energyZ = stackEnergy.get(str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j]), 1)
     
        # Calculate the Boltzmann factor and scale by 2, then multiply by the Zb of the next pair
        #print("stack multiple of", math.exp(-beta*(energyG - scale*scale_num)))
        if debug_energy: print("stack of {}, {} on {}, {}:".format(s[i], s[j], s[i+1], s[j-1]), (-1/beta)*math.log(energyZ))
        return energyZ*traceZbd(i+1)

    # Internal loop case
    if case == 1:
        d = get_next_loop(i)
        e = pair_dict[d]

        # Set the energy values
        energyZ = 1
                        
        #print("e is %d"%e)
        if e >= j:
            raise Exception("Went out of bounds. This is not good")
            return -1
                
        # Break if loop is too big
        if (d-i-1) + (j-e-1) >= l:
            raise Exception("Went out of bounds. This is not good")
            return -1

        # The scale number is (d-i) + (j-e)
        scale_num = (d-i) + (j-e)
        # Number in top loop
        top = d-i-1
        # Number in bottom loop
        bottom = j-e-1
                            
        # Base situation (we don't want to count this)
        if top == 0 and bottom == 0:
            raise Exception("this should have been handled by the stacking case")
            return -1
        # Bulge situation
        elif top == 0 or bottom == 0:
            energyZ *= bulgeLoopEnergy[top+bottom]
            #print("bulge energy is", bulgeLoopEnergy[top+bottom])
            # End Penalty
            if top > 1 or bottom > 1:
                #print(energyG)
                energyZ *= endPenalty.get(str(s[i])+str(s[j]), 1)
                #print("auPen1:", endPenalty.get(str(s[i])+str(s[j]), 0))
                energyZ *= endPenalty.get(str(s[e])+str(s[d]), 1)
                #print("auPen2:", endPenalty.get(str(s[e])+str(s[d]), 0)) 
            # Stacking energy (added 5 Dec 2016)
            if top == 1 or bottom == 1:
                energyZ *= stackEnergy.get(str(s[i])+str(s[d])+str(s[e])+str(s[j]), 1)
                #print(stackEnergy.get(str(s[i])+str(s[d])+str(s[e])+str(s[j]), 0))
            #print("traceback added", energyG)
        # 1x1 situation
        elif top == 1 and bottom == 1:
            energyZ *= int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])]
        # 2x1 situation
        elif top == 2 and bottom == 1:
            energyZ *= int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])]
        # 1x2 situation
        elif top == 1 and bottom == 2:
            energyZ *= int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
        # 2x2 situation
        elif top == 2 and bottom == 2:
            energyZ *= int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
        # General case
        # See http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
        else:
            #print("else")
            # Initiation term
            energyZ *= internalLoopEnergy[top+bottom]
            #print("size", top+bottom, "->", energyG)
        
            # Asymmetry term
            # Coef rounded from 0.48 to 0.5 to follow unafold and RNAeval with 99 rules
            energyZ *= asymmetry[abs(top-bottom)]
            #print("asymmetry", min(3, 0.5*abs(top-bottom)))
        
            # 'First Mismatch' term
            if top > 1 and bottom > 1:
                # 5' side
                energyZ *= loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 1)
                #print("5' 1st", loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 0))
                # 3' side
                energyZ *= loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 1)
                #print("3' 1st", loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 0))

            # Loop closure penalty term
            energyZ *= loopEndPenalty.get(str(s[i])+str(s[j]), 1)
            energyZ *= loopEndPenalty.get(str(s[e])+str(s[d]), 1)
            #print("5' AU", loopEndPenalty.get(str(s[i])+str(s[j]), 0))
            #print("3' AU", loopEndPenalty.get(str(s[e])+str(s[d]), 0))
                  
        if debug_energy: print("internal/bulge loop multiple of", energyZ)
        return energyZ*traceZbd(d)
                
    # Multibranch case
    # Remember, the dangles are backwards in this case
    if case < 2:
        raise Exception("This was supposed to be the only other case. Something is wrong")
    #print("case 2")
    ### WILL NEED TO INCLUDE THE SCALE FACTOR ###
    
    endPen = endPenalty.get(str(s[i])+str(s[j]), 1)
    # No dangle
    if not i in dangle_list and not j-1 in dangle_list:
        if debug_energy:
            print("multibranch ea:", (-1/beta)*math.log(ea))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch ec:", (-1/beta)*math.log(ec))
                    
        return (ea*ec
               * traceZ1d(i+1, j-1, case)
               * endPen)
    
    # 3'D
    if i in dangle_list and not j-1 in dangle_list:
        if debug_energy: 
            print("multibranch ea:", (-1/beta)*math.log(ea))
            print("multibranch eb:", (-1/beta)*math.log(eb))
            print("multibranch ec:", (-1/beta)*math.log(ec))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch 3dangle:", (-1/beta)*math.log(dangle3Energy[str(s[j])+str(s[i])+str(s[i+1])]))

        return (ea*eb*ec*dangle3Energy[str(s[j])+str(s[i])+str(s[i+1])]
               * traceZ1d(i+2, j-1, case)
               * endPen)
    
    # 5'D
    if not i in dangle_list and j-1 in dangle_list:
        if debug_energy: 
            print("multibranch ea:", (-1/beta)*math.log(ea))
            print("multibranch eb:", (-1/beta)*math.log(eb))
            print("multibranch ec:", (-1/beta)*math.log(ec))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch 5dangle:", (-1/beta)*math.log(dangle5Energy[str(s[j-1])+str(s[j])+str(s[i])]))

        return (ea*eb*ec*dangle5Energy[str(s[j-1])+str(s[j])+str(s[i])]
               * traceZ1d(i+1, j-2, case)
               * endPen)
    
    # ts
    if i in dangle_list and j-1 in dangle_list:
        #if (s[i+1], s[j-1]) in allowedPairs:
        #    print("they are {},{}".format(i+1, j-1))
        #    print("aka {},{}".format(s[i+1], s[j-1]))
        #    raise Exception("this was not supposed to happen")
            #return -1
        if debug_energy: 
            print("multibranch ea:", (-1/beta)*math.log(ea))
            print("multibranch eb*eb:", (-1/beta)*math.log(eb*eb))
            print("multibranch ec:", (-1/beta)*math.log(ec))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch ts:", (-1/beta)*math.log(terminalMismatchMEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]))

        return (ea*eb*eb*ec*terminalMismatchMEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]
               * traceZ1d(i+2, j-2, case)
               * endPen)

    raise Exception("Did not reach random number in probability sum")
    return -1 

# Traceback for Z1
def traceZ1d(i, j, num_structures):
    if debug_structure: print("traceZ1d({},{},{})".format(i, j, num_structures))
    # Base case
    if i == j: 
        raise Exception("something is wrong")
        return -1 
    
    # Error case
    if i > j:
        raise ValueError("i greater than j, something is wrong")
        return "Something went wrong"
    
    # We have a few possibilities.
    # i unpaired    
    # i, k paired, later structure
    # i, k paired, no later structure

    # i unpaired
    if not i in pair_dict and not (i+1 in pair_dict and i in dangle_list):
        if debug_energy: print("multibranch eb:", (-1/beta)*math.log(eb))
        return eb*traceZ1d(i+1, j, num_structures)
    
    # i, k paired, then no more structure
    if num_structures == 1:
        # i, k paired
        if i in pair_dict:
            k = pair_dict[i]
            endPen = endPenalty.get(str(s[i])+str(s[k]), 1)

            # no dangle
            if not i-1 in dangle_list and not k in dangle_list:
                scale_num = j-k
                if debug_energy:
                    print("multibranch ec:", (-1/beta)*math.log(ec))
                    print("AU penalty:", (-1/beta)*math.log(endPen))
                    print("multibranch eb**(j-k):", (-1/beta)*math.log(eb**(j-k)))
                return (ec
                       * traceZbd(i)
                       * (eb**(j-k))
                       * endPen)

            # 3'D
            if not i-1 in dangle_list and k in dangle_list:
                if k > min(j-1, N-2): 
                    raise Exception("oof")
                scale_num = j-k
                if debug_energy:
                    print("multibranch ec:", (-1/beta)*math.log(ec))
                    print("AU penalty:", (-1/beta)*math.log(endPen))
                    print("multibranch eb**(j-k):", (-1/beta)*math.log(eb**(j-k)))
                    print("multibranch 3dangle:", (-1/beta)*math.log(dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]))
                return (ec
                       * traceZbd(i)
                       * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                       * (eb**(j-k))
                       * endPen)

        # i+1, k paired
        else:
            k = pair_dict[i+1]
            endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)

            # 5'D
            if i in dangle_list and not k in dangle_list:
                if k > min(j, N-1):
                    raise Exception("oof")
                scale_num = j-k+1
                if debug_energy: 
                    print("multibranch ec:", (-1/beta)*math.log(ec))
                    print("AU penalty:", (-1/beta)*math.log(endPen))
                    print("multibranch eb**(j-k+1):", (-1/beta)*math.log(eb**(j-k+1)))
                    print("multibranch 5dangle:", (-1/beta)*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]))
                return (ec
                       * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
                       * traceZbd(i+1)
                       * (eb**(j-k+1))
                       * endPen)
            
            # ts
            if i in dangle_list and k in dangle_list:
                if k > min(j-1, N-2):
                    raise Exception("nope")
                scale_num = j-k+1
                # Here we treat the double dangle energy as a external terminal mismatch energy ('tstackm')
                # ^^ We will check if this is correct
                #if (s[i], s[k+1]) in allowedPairs:
                    #raise Exception("what")
                if debug_energy: 
                    print("multibranch ec:", (-1/beta)*math.log(ec))
                    print("AU penalty:", (-1/beta)*math.log(endPen))
                    print("multibranch eb**(j-k+1):", (-1/beta)*math.log(eb**(j-k+1)))
                    print("multibranch ts:", (-1/beta)*math.log(terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]))
                return (ec
                       * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
                       * traceZbd(i+1)
                       * (eb**(j-k+1))
                       * endPen)

    # i, k paired, some later structure
    if num_structures < 2:
        raise Exception("this case should have been handled")
    #print("num_structures is big", num_structures)

    if i in pair_dict:
        k = pair_dict[i]
        endPen = endPenalty.get(str(s[i])+str(s[k]), 1)
        #print("k is", k)

        # No dangles
        if not i-1 in dangle_list and not k in dangle_list:
            if debug_energy:
                print("multibranch ec:",(-1/beta)*math.log(ec))
            return (ec
                   * traceZbd(i)
                   * traceZ1d(k+1, j, num_structures - 1)
                   * endPen)
    
        # 3'D
        if not i-1 in dangle_list and k in dangle_list:
            if k > min(j-2, N-3):
                raise Exception("th")
            if debug_energy: 
                print("multibranch eb", (-1/beta)*math.log(eb))
                print("multibranch ec:", (-1/beta)*math.log(ec))
                print("AU penalty:", (-1/beta)*math.log(endPen))
                print("multibranch 3dangle:", (-1/beta)*math.log(dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]))
            return (eb*ec
                   * traceZbd(i)
                   * dangle3Energy[str(s[i])+str(s[k])+str(s[k+1])]
                   * traceZ1d(k+2, j, num_structures - 1)
                   * endPen)
    
    # i+1, k paired
    k = pair_dict[i+1]
    endPen = endPenalty.get(str(s[i+1])+str(s[k]), 1)
    #if k == 19: print("ok1")
    #if (i+1, k) in pair_dict.items(): print("ok2")

    # 5'D
    if i in dangle_list and not k in dangle_list:
        if k > min(j-1, N-2):
            raise Exception("lu")
        if debug_energy: 
            print("multibranch eb", (-1/beta)*math.log(eb))
            print("multibranch ec:", (-1/beta)*math.log(ec))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch 5dangle:", (-1/beta)*math.log(dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]))
        return (eb*ec
               * traceZbd(i+1)
               * dangle5Energy[str(s[i])+str(s[i+1])+str(s[k])]
               * traceZ1d(k+1, j, num_structures - 1)
               * endPen)
    
    # ts
    if i in dangle_list and k in dangle_list:
        if k > min(j-2, N-3): 
            raise Exception("lu2")
        #if (s[i], s[k+1]) in allowedPairs:
            #raise Exception("bad")
        if debug_energy:
            print("multibranch eb*eb", (-1/beta)*math.log(eb*eb))
            print("multibranch ec:", (-1/beta)*math.log(ec))
            print("AU penalty:", (-1/beta)*math.log(endPen))
            print("multibranch ts:", (-1/beta)*math.log(terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]))
        return (eb*eb*ec
               * traceZbd(i+1)
               * terminalMismatchMEnergy[str(s[k])+str(s[k+1])+str(s[i])+str(s[i+1])]
               * traceZ1d(k+2, j, num_structures - 1)
               * endPen)
    
    # If we get here, then something is wrong with our probabilities, 
    raise Exception("Did not reach random number in probability sum")
    return -1


##################
#  Zb functions  #
##################

################
# Forward Fill #
################
# The following functions are only ever called during the forward fill. There are separate
# functions to be called during the wraparound.

# Get the hairpin term
# i is paired with j, the hairpin starts immediately after
def hairpinTerm(i, j):
    #print("Pair {},{}:".format(i, j))

    # Calculated according to http://rna.urmc.rochester.edu/NNDB/turner99/hairpin.html
    # How many bases in hairpin (not including final stack)
    num = j - i - 1
    # How many bases to count for scale factor
    scale_num = num + 2
    # The sequence in the hairpin (from i to j, inclusive)
    seq = s[i: j+1]
    #print(seq)
    
    # The energy to return
    energyZ = 1
    
    # If there are not enough bases to make a hairpin, return 0 energy
    # We also consider it impossible to have more than l (usually 30) bases in our hairpin
    if num < 3: return energyZ
    
    # There is a separate (simplified) procedure when there are only three bases in the hairpin
    elif num == 3:
        # We set the total energy to the initiation energy of the hairpin
        # Only other factor to consider is C-loops, which are added on after the if/else
        energyZ = hairpinEnergy[num]
        #print("Size: 3, Energy:", energy, "->", math.exp(-beta*energy))
        # Add on the AU penalty
        energyZ *= endPenalty.get(str(s[i])+str(s[j]), 1)
        #if endPenalty.get(str(s[i])+str(s[j]), 0) > 0: print("    with endPen:", energy, "->", math.exp(-beta*energy))
        #print("AU:", endPenalty.get(str(s[i])+str(s[j]), 0))

    # The 'usual' case, where there are at least four bases in the hairpin
    else:
        #ΔG°37 total =      
        #ΔG°37 initiation (n) +
        #ΔG°37 (terminal mismatch) +
        #ΔG°37 (UU or GA first mismatch) +
        #ΔG°37 (special GU closure) +
        #ΔG°37 penalty (all C loops)

        # NOTE: following along to the UNAfold output,
        # we add in the endPenalty of 0.5 like in other places
        # (done on 4 Dec)

        # NOTE: the website lists terminal mismatch and first mismatch
        # as two separate entities to be added idividually in the energy sum.
        # As done in UNAfold, this MacroFold implementation
        # uses terminalMismatchHEnergy (tstackh.csv) values
        # to account for both of these energy penalties/bonuses
        # simultaneously. This is why neither of those
        # terms appear in the sum below.
        # (noted on 4 Jan when testing agaist RNAsubopt)
    
        # We set the total energy to the initiation energy of the hairpin,
        # and then consider and add all other factors
        if num <= 30:
            energyZ = hairpinEnergy[num]
        else:
            energyZ = math.exp(-beta*(6.4 + 1.75*r*T*math.log(num/9)))
        #print("Size: {}, Energy:".format(num), energy, "->", math.exp(-beta*energy))
        #print("Size:", energy)
        #print("energy is", energy)
    
        # Add on the terminal mismatch energy
        energyZ *= terminalMismatchHEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]
        #if terminalMismatchHEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])] > 0:
            #print("    with endPen:", energy, "->", math.exp(-beta*energy))

    # Add on the special GU closure energy
    # Check we don't fall off the front
    if i-2 >= 0:
        energyZ *= guClosureEnergy.get(str(s[i-2])+str(s[i-1])+str(s[i])+str(s[j]), 1)

    # Add on special loop bonuses
    energyZ *= specialLoops.get(seq, 1)
    #print(specialLoops)
    #if seq in specialLoops:
        #print("WAHOO, {} is a tasty special loop!".format(seq))
        
    # Add on the C-loop bonus (only considers the internal parts)
    energyZ *= CLoops.get(seq[1:-1], 1)
    
    #print("Total Hairpin Energy:  %s"%energy)
    #print("Initiation:        %s"%hairpinEnergy[num])
    #if num > 3:
        #if terminalMismatchHEnergy.get(str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j]), 0) != 0:
            #print("Terminal Mismatch: %s"%terminalMismatchHEnergy.get(str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j]), 0))
        #if specialLoops.get(seq, 0) != 0:
            #print("Special Loops:     %s"%specialLoops.get(seq, 0))
        #if CLoops.get(seq[1:-1], 0) != 0:
            #print("C-loops:           %s"%CLoops.get(seq[1:-1], 0))
    #print("total:", energy)
    #print("-> ", math.exp(-beta*energy))
        
    #print()
    
    # Calculate the Boltzmann factor and scale by number of included bases
    return energyZ

# Get the stack term
# i is paired with j (already checked to make sure it is legal)
def stackTermF(i, j):
    
    # The scale number is 2 because it's a stack
    scale_num = 2
    
    # Calculated according to http://rna.urmc.rochester.edu/NNDB/turner99/wc.html
    # Calculate the stacking energy of i, j with i+1, j-1
    # We do not need to worry about the case in which it cannot stack, as that will be taken
    # care of by the fact that in the end we multiply by Zb[i+1, j-1]
    energyZ = stackEnergy.get(str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j]), 1)
    #if Zb[i+1, j-1] != 0:
        #print("total:", energy)
        #print("->", math.exp(-beta*energy))
        #print("to be multiplied by", Zb[i+1, j-1]) 
        #print("->", math.exp(-beta*(energy - scale*scale_num))*Zb[i+1, j-1])
 
    # Printout for testing purposes
    #if energy != 0:
    #    print("stack energy for %s%s"%(s[i], s[i+1]))
    #    print("                 %s%s is: %s"%(s[j], s[j-1], energy))
    #    print("will be multiplied by %s"%Zb[i+1, j-1])
    
    # Calculate the Boltzmann factor and scale by 2, then multiply by the Zb of the next pair
    return energyZ*Zb[i+1, j-1]
    
# Get the multibranch term
# i is paired with j (already checked to make sure it is legal)
def multibranchTermF(i, j):
    
    endPen = endPenalty.get(str(s[i])+str(s[j]), 1)
    # The scale number is 2 because it's a multibranch
    scale_num = 2
    
    # Set the return energy
    energyZ = 0
    
    # No dangles
    energyZ += ea*ec*Z2[i+1, j-1]*endPen
    
    # 3'D
    energyZ += ea*eb*ec*dangle3Energy[str(s[j])+str(s[i])+str(s[i+1])]*Z2[i+2, j-1]*endPen
    
    # 5'D
    energyZ += ea*eb*ec*dangle5Energy[str(s[j-1])+str(s[j])+str(s[i])]*Z2[i+1, j-2]*endPen
    
    # ts
    #if 1:#(s[i+1], s[j-1]) not in allowedPairs:
    if not (s[i+1], s[j-1]) in allowedPairs:
        # ^ why did we do this?
        # I don't know but I switched it back (16 Jan)
        energyZ += ea*eb*eb*ec*terminalMismatchMEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]*Z2[i+2, j-2]*endPen
    
    # Calculate the Boltzmann factor and scale by 2, then multiply by the Z2 of the next pair
    return energyZ

# Get the internal loop term
# i is paired with j
def internalTermF(i, j):
    #print("looking at %s,%s pair"%(i, j))
    # What is our allowed internal loop size
    global l
    global e_index
    
    # Border case
    if j-i < 7: return 0 # changed from <= to < on 11 Jan
    
    # The implementation of internal loop calculations is rather complex
    # We follow that found at http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
    
    # Here we loop over the d's less than l away from i, and then over the e's less than l away from j
    #for dis in range(4, j-i-1):
    #    d = i+1
    #    while d in range(i+1, i+l+1) and d < N-5: # need to also check to make sure we're not past j-1
    #        e = d + dis
    #for i in range(N-4):
    #    for delta in range(5,N-4-i):
    #        j=i+delta

    # The total energy, in partition function form
    energyZ_tot = 0

    for d in range(i+1, min(i+l,j-4)): # changed from j to j-4 (Jan 11)
        # We use this variable for the duration of this loop
        if not (d, j) in e_index: continue

        eIdx = e_index[(d, j)]
        # Get the rightmost possible e
        e = R[d][eIdx]
        #print("d is", d)
        #print("e is", e)
        #print("j is", j)
        
        while (j-e-1) + (d-i-1) < l and eIdx >= 0:

            # We don't allow it when the end breaks the loop structure
            #if (i < N and j >= N) and (d >= N or e < N): 
                #return 0
            # ^ I don't think we ever see this in the forward-fill? (Jan 11)

            #print("eIdx is ", eIdx)
            #print("d, e are {}, {}".format(d, e))
            #print("we are at d = {}".format(d))
            #print("R[{}] = {}".format(d, R[d]))
            
            # The energy of this particular configuration
            energyZ = 1

            #print("e is %d"%e)
            if e >= j: break
                
            # Break if loop is too big
            if d-i-1 + j-e-1 >= l: break # as of 20 Jan this should never be called
            
            # The scale number is (d-i) + (j-e)
            scale_num = (d-i) + (j-e)
            # Number in top loop
            top = d-i-1
            # Number in bottom loop
            bottom = j-e-1

            #print(olds)
            #print(' '*i+'^'+' '*(d-i-1)+'*'+' '*(e-d-1)+'*'+' '*(j-e-1)+'^')
            #print("R[d] is", R[d])
            #print("starting eIdx is", eIdx)
            #print("e is", e)

            # Base situation (we don't want to count this)
            if top == 0 and bottom == 0:
                eIdx -= 1
                e = R[d][eIdx]
                #print("added no energy")
                continue
            # Bulge situation
            elif top == 0 or bottom == 0:
                energyZ *= bulgeLoopEnergy[top+bottom]
                #print("bulge energy is", bulgeLoopEnergy[top+bottom])
                # End Penalty
                if top > 1 or bottom > 1:
                    energyZ *= endPenalty.get(str(s[i])+str(s[j]), 1)
                    energyZ *= endPenalty.get(str(s[e])+str(s[d]), 1)
                # Stacking energy (added 5 Dec 2016, "if" added 11 Jan 2017)
                if top == 1 or bottom == 1:
                    energyZ *= stackEnergy.get(str(s[i])+str(s[d])+str(s[e])+str(s[j]), 1)
                #print("forwardfill added", energyG)
            # 1x1 situation
            elif top == 1 and bottom == 1:
                energyZ *= int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])]
                #print("int11 for ", str(s[i])+str(s[i+1])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[e]))
                #print("is: ", int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])])
            # 2x1 situation
            elif top == 2 and bottom == 1:
                energyZ *= int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])]
                #print("int21 for ", str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[e]))
                #print("is: ", int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])])
            # 1x2 situation
            elif top == 1 and bottom == 2:
                energyZ *= int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
                #print("int21 for ", str(s[i])+str(s[i+1])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[j-2])+str(s[e]))
                #print("is: ", int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])])
            # 2x2 situation
            elif top == 2 and bottom == 2:
                energyZ *= int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
                #print("int22 for ", str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[j-2])+str(s[e]))
                #print("is: ", int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])])
            # General case
            # See http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
            else:
                #print("Loop: {}".format(s[i:d+1]))
                #print("      {}".format(s[e:j+1][::-1]))
                # Initiation term
                energyZ *= internalLoopEnergy[top+bottom]
                #print("Internal Energy: %s"%internalLoopEnergy[top+bottom])
                
                # Asymmetry term
                # Coef rounded from 0.48 to 0.5 to follow unafold and RNAeval with 99 rules
                energyZ *= asymmetry[abs(top-bottom)]
                #print("Asymmetry: %s"%(0.5*abs(top-bottom)))
                
                # 'First Mismatch' term
                if top > 1 and bottom > 1:
                    # 5' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 1)
                    # 3' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 1)
                    #print("mismatches: %s, %s"%(loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 0),
                    #                            loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 0)))
                
                # Loop closure penalty term
                energyZ *= loopEndPenalty.get(str(s[i])+str(s[j]), 1)
                energyZ *= loopEndPenalty.get(str(s[e])+str(s[d]), 1)
                #print("closures: %s, %s"%(loopEndPenalty.get(str(s[i])+str(s[j]), 0),
                #                          loopEndPenalty.get(str(s[e])+str(s[d]), 0)))
                
            # Add this case to the total energy in exponential form
            energyZ_tot += energyZ*Zb[d, e]
            #if energyZ > 1:
                #print("aha1")
                #print(olds)
                #raise Exception("{},{},{},{}".format(i, d, e, j))
                    
            # Change the e_index
            eIdx -= 1
            #print("new ")

            # We get the actual value of e
            e = R[d][eIdx]

        #if energyZ > 1:
            #print("aha2")
    #if energyZ > 1:
        #print("aaagh")
        #raise Exception("{},{},{},{}".format(i, d, e, j))
    return energyZ_tot

########################
# Wraparound Functions #
########################
# These are only called in the wraparound case. Thus, for these functions:
# * i < N <= j < min(2N, i+N-1)
# * s = 2*s


# Get the external loop term
# i is paired with j
# QUESTION: is there really no au penalty term
def externalLoopTerm(i, j):

    if i == N-1 and j == N:
        #print("filled bottom left")
        return 1
    elif i == N-1 and j == N+1:
        return 1 + dangle5Energy[str(s[0])+str(s[1])+str(s[N-1])]
    elif i == N-1:
        return Z5[j-N-1] + Z5[j-N-2]*dangle5Energy[str(s[j-N-1])+str(s[j-N])+str(s[N-1])]
    elif i == N-2 and j == N:
        return 1 + dangle3Energy[str(s[0])+str(s[i])+str(s[i+1])]
    elif j == N:
        return Z3[i+1] + Z3[i+2]*dangle3Energy[str(s[0])+str(s[i])+str(s[i+1])]
    elif i == N-2 and j == N+1:
        energyZ = 1 + dangle3Energy[str(s[j-N])+str(s[i])+str(s[i+1])] + dangle5Energy[str(s[j-N-1])+str(s[j-N])+str(s[i])]
        if not (s[i+1], s[j-N-1]) in allowedPairs:
            energyZ += terminalMismatchEEnergy[str(s[i])+str(s[i+1])+str(s[j-N-1])+str(s[j-N])]
        return energyZ
    elif i == N-2:
        energyZ = (Z5[j-N-1]
                  + Z5[j-N-1]*dangle3Energy[str(s[j-N])+str(s[i])+str(s[i+1])]
                  + Z5[j-N-2]*dangle5Energy[str(s[j-N-1])+str(s[j-N])+str(s[i])])
        if not (s[i+1], s[j-N-1]) in allowedPairs:
            energyZ += Z5[j-N-2]*terminalMismatchEEnergy[str(s[i])+str(s[i+1])+str(s[j-N-1])+str(s[j-N])]
        return energyZ
    elif j == N+1:
        energyZ = (Z3[i+1]
                  + Z3[i+2]*dangle3Energy[str(s[j-N])+str(s[i])+str(s[i+1])]
                  + Z3[i+1]*dangle5Energy[str(s[j-N-1])+str(s[j-N])+str(s[i])])
        if not (s[i+1], s[j-N-1]) in allowedPairs:
            energyZ += Z3[i+2]*terminalMismatchEEnergy[str(s[i])+str(s[i+1])+str(s[j-N-1])+str(s[j-N])]
        return energyZ
    else:
        # consider the four cases: no dangle, 3'D, 5'D, and tstack
        # No dangle
        energyZ = Z3[i+1]*Z5[j-N-1]
            
        # 3'D
        #print("i is", i)
        #print("i+2 is", i+2)
        energyZ += Z3[i+2]*Z5[j-N-1]*dangle3Energy[str(s[j-N])+str(s[i])+str(s[i+1])]
        
        # 5'D
        energyZ += Z3[i+1]*Z5[j-N-2]*dangle5Energy[str(s[j-N-1])+str(s[j-N])+str(s[i])]
        
        # tstack
        if not (s[i+1], s[j-N-1]) in allowedPairs:
            energyZ += Z3[i+2]*Z5[j-N-2]*terminalMismatchEEnergy[str(s[i])+str(s[i+1])+str(s[j-N-1])+str(s[j-N])]

        return energyZ

# Get the stack term
# i is paired with j (already checked to make sure it is legal)
def stackTermW(i, j):

    # Base case: we're too close to the ends to stack
    if i == N-1 or j == N:
        return 0

    # The scale number is 2 because it's a stack
    scale_num = 2
    
    # Calculated according to http://rna.urmc.rochester.edu/NNDB/turner99/wc.html
    # Calculate the stacking energy of i, j with i+1, j-1
    # We do not need to worry about the case in which it cannot stack, as that will be taken
    # care of by the fact that in the end we multiply by Zb[i+1, j-1]
    energy = stackEnergy.get(str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j]), 1)
 
    # Printout for testing purposes
    #if energy != 0:
    #    print("stack energy for %s%s"%(s[i], s[i+1]))
    #    print("                 %s%s is: %s"%(s[j], s[j-1], energy))
    #    print("will be multiplied by %s"%Zb[i+1, j-1])
    
    # Calculate the Boltzmann factor and scale by 2, then multiply by the Zb of the next pair
    #if i == 5 and j == 14: print("stack", math.exp(-beta*(energy - scale*scale_num))*Zb[i+1, j-1-N])
    return energy*Zb[i+1, j-1-N]

# Get the multibranch term
# i is paired with j (already checked to make sure it is legal)
def multibranchTermW(i, j):
    #if i == 6 and j == 15:
        #print("yeehaw")
    endPen = endPenalty.get(str(s[i])+str(s[j]), 1)

    # Base cases
    if i >= N-1 or j <= N:
    #if i >= N or j < N: # changed on 18 Jan
        #print(i, j, "hmm")
        return 0

    # Set the j to where it would have been 
    j = j-N

    # Then from here on it's exactly the same as the usual forward fill
    
    # The scale number is 2 because it's a multibranch
    scale_num = 2
    
    # Set the return energy
    energyZ = 0
    
    # No dangles
    energyZ += Z2[i+1, j-1]*endPen
    #if i == 6 and j == 15-N: print("no dang energyZ", energyZ)

    # 3'D
    ## CHECK DIRECTION OF THIS AND NEXT DANGLE ##
    if i <= N-3:
        energyZ += eb*dangle3Energy[str(s[j])+str(s[i])+str(s[i+1])]*Z2[i+2, j-1]*endPen
    #if i == 6 and j == 15-N: print("3 dang energyZ", energyZ)
    
    # 5'D
    #if j >= N+2:
    if j >= 2: # changed on 19 Jan
        energyZ += eb*dangle5Energy[str(s[j-1])+str(s[j])+str(s[i])]*Z2[i+1, j-2]*endPen
    #if i == 6 and j == 15-N: print("5 dang energyZ", energyZ)
    
    # ts
    #if 1:#(s[i+1], s[j-1]) not in allowedPairs:
    #if i < N-2 and j > N+2: changed on 19 Jan
    if i <= N-3 and j >= 2:
        if not (s[i+1], s[j-1]) in allowedPairs:
            energyZ += eb*eb*terminalMismatchMEnergy[str(s[i])+str(s[i+1])+str(s[j-1])+str(s[j])]*Z2[i+2, j-2]*endPen
            #if energyZ != 0: print("energyZ", energyZ)
    
    # Calculate the Boltzmann factor and scale by 2, then multiply by the Z2 of the next pair
    #if i == 0 and j == 11:
    #if i == 6 and j == 15-N: print("returning", math.exp(-beta*(ea + ec))*energyZ)
    #if i == 6 and j == 15-N: print("we see", Z2[7,14-N])
    return ea*ec*energyZ

# Get the internal loop term
# i is paired with j
def internalTermW(i, j):

    # Base case
    if i >= N-1 or j <= N or j >= i+N-3: #(replaced i >= N-2 with i >= N-1) Jan 18
        #print("edge")
        return 0

    #print("looking at %s,%s pair"%(i, j))
    # What is our allowed internal loop size
    global l
    global e_index
    
    # The implementation of internal loop calculations is rather complex
    # We follow that found at http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
    
    # Case energy in G and Z forms
    energyZ_tot = 0
    
    # Here we loop over the d's less than l away from i, and then over the e's less than l away from j
    #for dis in range(4, j-i-1):
    #    d = i+1
    #    while d in range(i+1, i+l+1) and d < N-5: # need to also check to make sure we're not past j-1
    #        e = d + dis
    #for i in range(N-4):
    #    for delta in range(5,N-4-i):
    #        j=i+delta
            
    for d in range(i+1, min(i+l,j)):
        #print(d)
                    
        # Stop when we reach the break point
        if d >= N: 
            break

        # We use this variable for the duration of this loop
        if not (d, j) in e_index: 
            #print("{},{} not in e_index".format(d,j))
            continue
        #print("d is ", d)
        #print("R[d] is ", R[d])
        eIdx = e_index[(d, j)]
        #print("eIdx", eIdx)
        #print("starting eIdx is ", eIdx)
        e = R[d][eIdx]
        #print("starting e is ", e)
        
        while (j-e-1) + (d-i-1) < l and eIdx >= 0 and e >= N:

            #print("eIdx is ", eIdx)
            #print("d, e are {}, {}".format(d, e))
            #print("we are at d = {}".format(d))
            #print("R[{}] = {}".format(d, R[d]))
            
            # We get the actual value of e
            #print(i, d, e, j)

            # Reset the energy of this microstate
            energyZ = 1

            #print("e is %d"%e)
            if e >= j: break
                
            # Break if loop is too big
            if d-i-1 + j-e-1 >= l: break

            # The scale number is (d-i) + (j-e)
            scale_num = (d-i) + (j-e)
            # Number in top loop
            top = d-i-1
            # Number in bottom loop
            bottom = j-e-1
                
            #print(i, d, e, j)
            # Base situation (we don't want to count this)
            if top == 0 and bottom == 0:
                #print("stack")
                eIdx -= 1
                e = R[d][eIdx]
                continue
            # Bulge situation
            elif top == 0 or bottom == 0:
                energyZ *= bulgeLoopEnergy[top+bottom]
                # End Penalty
                if top > 1 or bottom > 1:
                    energyZ *= endPenalty.get(str(s[i])+str(s[j]), 1)
                    energyZ *= endPenalty.get(str(s[e])+str(s[d]), 1)
                # Stacking energy (added 5 Dec 2016)
                if top == 1 or bottom == 1:
                    energyZ *= stackEnergy.get(str(s[i])+str(s[d])+str(s[e])+str(s[j]), 1)
            # 1x1 situation
            elif top == 1 and bottom == 1:
                energyZ *= int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])]
                #print("int11 for ", str(s[i])+str(s[i+1])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[e]))
                #print("is: ", int11[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-1])+str(s[j])])
            # 2x1 situation
            elif top == 2 and bottom == 1:
                energyZ *= int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])]
                #print("int21 for ", str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[e]))
                #print("is: ", int21[str(s[e])+str(s[j-1])+str(s[j])+str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])])
            # 1x2 situation
            elif top == 1 and bottom == 2:
                energyZ *= int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
                #print("int21 for ", str(s[i])+str(s[i+1])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[j-2])+str(s[e]))
                #print("is: ", int21[str(s[i])+str(s[i+1])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])])
            # 2x2 situation
            elif top == 2 and bottom == 2:
                energyZ *= int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])]
                #print("int22 for ", str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d]))
                #print("on        ", str(s[j])+str(s[j-1])+str(s[j-2])+str(s[e]))
                #print("is: ", int22[str(s[i])+str(s[i+1])+str(s[i+2])+str(s[d])+str(s[e])+str(s[j-2])+str(s[j-1])+str(s[j])])
            # General case
            # See http://rna.urmc.rochester.edu/NNDB/turner99/internal.html
            else:
                #print("Loop: {}".format(s[i:d+1]))
                #print("      {}".format(s[e:j+1][::-1]))
                # Initiation term
                energyZ *= internalLoopEnergy[top+bottom]
                #print("Internal Energy: %s"%internalLoopEnergy[top+bottom])
                
                # Asymmetry term
                # Coef rounded from 0.48 to 0.5 to follow unafold and RNAeval with 99 rules
                energyZ *= asymmetry[abs(top-bottom)]
                #print("Asymmetry: %s"%(0.5*abs(top-bottom)))
                
                # 'First Mismatch' term
                if top > 1 and bottom > 1:
                    # 5' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 1)
                    # 3' side
                    energyZ *= loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 1)
                    #print("mismatches: %s, %s"%(loopFirstMismatchEnergy.get(str(s[i+1])+str(s[j-1]), 0),
                    #                            loopFirstMismatchEnergy.get(str(s[e+1])+str(s[d-1]), 0)))
                
                # Loop closure penalty term
                energyZ *= loopEndPenalty.get(str(s[i])+str(s[j]), 1)
                energyZ *= loopEndPenalty.get(str(s[e])+str(s[d]), 1)
                #print("closures: %s, %s"%(loopEndPenalty.get(str(s[i])+str(s[j]), 0),
                #                          loopEndPenalty.get(str(s[e])+str(s[d]), 0)))
                
            # Add this case to the total energy in exponential form
            energyZ_tot += energyZ*Zb[d, e-N]
            #if i == 5 and j == 14:
                #print("internal", math.exp(-beta*(energy - scale*scale_num))*Zb[d, e-N])
                    
            # Change the e_index
            eIdx -= 1

            # Get the new e
            e = R[d][eIdx]
            #print("new ")
    return energyZ_tot


###########
# Testing #
###########

# Enter debug mode
def debug_mode():
    global debug_structure
    global debug_energy

    debug_structure = True
    debug_energy = True

# Generate a random sequence of length n
def generate_sequence(n):
    bases = ['A', 'C', 'G', 'U']
    seq = np.random.choice(bases, size=n)
    return(''.join(seq))

# have some kind of cache dict that we can store things in on the fly?

# Read the .ct into pair_dict and dangle_list
def read_ct(file_name):
    global pair_dict
    global dangle_list

    pair_dict = {}
    dangle_list = []

    # Open file as lines
    with open(file_name, 'r') as myfile:
        lines = myfile.readlines()

    # Split each line on whitespace
    for (i,line) in enumerate(lines):
        lines[i] = lines[i].split()
    # Remove the header
    lines = lines[1:]

    # Key columns are 4 and 6
    for (i,line) in enumerate(lines):

        # Add pair to the dict
        num = int(line[4]) 
        if num != 0:
            if i < num:
                pair_dict[i] = num - 1

        # Add non-stacking dangle
        if int(line[6]) != 0 and int(line[7]) == 0 and int(line[4]) == 0:
            dangle_list.append(i-1) 
        if int(line[7]) != 0 and int(line[6]) == 0 and int(line[4]) == 0:
            if not i in dangle_list: dangle_list.append(i)
    
    #print("pairs", pair_dict)
    #print("dangles", dangle_list)
    #print(lines)

# A general testing function
def run_tests(number=1000, length=40, hairpin=0):
    global debug_energy
    #debug_energy = True
    random = hairpin == 0
    #print(random)

    for i in range(number):
        if random: hairpin = generate_sequence(length)
        #print(hairpin)

        # Call the MacroFold function and print out energy of ends pairing
        loadSequence(hairpin)
        fillEnergyTables(no_dangles=False, zero_energies=False)
        fillAllowedPairsList()
        filleIndex()
        fillMatrices1()
        fillZ3()
        fillZ5()
        fillMatrices2()

        #print("          000000000011111111112222222222333333333344444444445555555555")
        #print("          012345678901234567890123456789012345678901234567890123456789")
        #print("Sequence: %s"%hairpin[:N])
        # I do believe this is the correct way to fill in the scale 
        #print("length is", N)
        # Print out the results
        #print("MacroFold: Zb = %s"%Zb[0][-1])
        #print("MacroFold: Z  = %s"%Z3[0])
        #with open('./test.dg') as file:
            #firstrow = True
            #for row in file:
                #if firstrow:
                    #firstrow = False
                    #continue
                #print("UNAFold:   Z = {}".format(row.split()[2]))
                #break
        #print()

        # Call the UNAFold function and print out same energy
        with open('./test.seq', 'w') as myfile:
            myfile.write(hairpin+'\n')
        remove = subprocess.getoutput("rm test.37.ct")
        subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/old/unafold-3.9/src/hybrid-ss --suffix=DAT -k 1 ./test.seq")
        output = subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/old/unafold-3.9/src/ct-energy --suffix=DAT ./test.37.ct") #-v -v ./test.37.ct")

        # Read in the .ct and assign the pair and dangle lists
        read_ct("test.37.ct")

        # Run the traceback and print the energy
        #print("UNAfold generated fold energy:", output)
        #print(subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/old/unafold-3.9/src/ct-energy --suffix=DAT -v -v ./test.37.ct"))
        #print("traceZ3d(0) is:", traceZ3d(0))
        output2 = -1/beta * math.log(traceZ3d(0))
        #print("%.4g" % output2)

        if abs(round(float(output), 1) - round(float(output2), 1)) < 1e-8:
            print(i, "agrees on", output, round(output2, 1))

        #elif abs(round(float(output), 1) - round(float(output2), 1)) - 0.5 < 1e-8:
            #print(i, "differs by 0.5")
        else:
            print("disagree:", output, "vs", output2)
            print("          000000000011111111112222222222333333333344444444445555555555")
            print("          012345678901234567890123456789012345678901234567890123456789")
            print("Sequence: %s"%hairpin[:N])
            debug_mode()
            print(subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/old/unafold-3.9/src/ct-energy --suffix=DAT -v -v ./test.37.ct"))
            traceZ3d(0)
            raise Exception("disagreed")
        #break

# Compare the partition functions
def compare_Z(sequence=False, length=14):

    if not sequence: sequence = generate_sequence(length)

    # Call the MacroFold function and print out energy of ends pairing
    loadSequence(sequence)
    fillEnergyTables()
    fillAllowedPairsList()
    filleIndex()
    fillMatrices1()
    fillZ3()
    fillZ5()
    fillMatrices2()

    print("          000000000011111111112222222222333333333344444444445555555555")
    print("          012345678901234567890123456789012345678901234567890123456789")
    print("Sequence: %s"%sequence[:N])
    #print("length is", N)
    '''print out the results'''
    #print("MacroFold: Zb = %s"%Zb[0][-1])
    #print("Zb =", Zb)
    print("MacroFold: Z = %s"%(-1/beta*math.log(Z3[0])))
    #print("Z3 =", Z3)
    with open('./test.seq', 'w') as myfile:
        myfile.write(sequence+'\n')
    remove = subprocess.getoutput("rm test.37.ct")
    #subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/unafold-3.9/src/hybrid-ss --suffix=DAT -k 1 ./test.seq")
    subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/unafold-3.9/src/hybrid-ss --suffix=DAT --prefilter=1 -E -t 37 -T 37 test.seq")
    subprocess.getoutput("RNAfold -v -p -d2 -P /usr/local/share/ViennaRNA/rna_turner1999.par < test.seq > seq1.out")
    # UNAfold
    with open('./test.dg') as myfile:
        firstrow = True
        for row in myfile:
            if firstrow:
                firstrow = False
                continue
            zed = float(row.split()[2])
            if zed > 1e-10:
                print("UNAFold:   Z = {}".format((-1/beta*math.log(float(row.split()[2])))))
            else:
                print("UNAFold:   Z = infinity")
            break
    # RNAfold
    with open("./seq1.out") as myfile:
        rows = myfile.read()
        #print(rows)
        for row in rows.split('\n'):
            if row.find('[') > -1:
                print("RNAfold:   Z = {}".format(row[row.find("[")+2:-1]))

# Get the individual states and calculate their energies
# with UNAfold, RNAfold, and MacroFold
def get_microstates(seq, zero_energies=False):
    global pair_dict
    global dangle_list

    seq = seq.upper()
    lock_n_load(seq, zero_energies=False)
    
    # Generate all the folds and put them in the folds folder
    #num_files = int(subprocess.getoutput("python3 ../scripts/get_all_structures.py {}".format(sequence)))
        
    num_folds, all_folds, all_structs = get_structs(str(seq), write=True)
    #for fold in all_folds:
        #print(fold)

    #print("num_folds is:", num_folds)

    # For each fold:
    energies = []
    differences = []
    total_Z = 0
    for i in range(num_folds):
        #print()
        #print("fold{}:".format(i+1))
        #print(seq)
        #print(all_folds[i])
        # Read in the .ct and assign the pair and dangle lists
        read_ct("/Users/niki/Google Drive/Thesis_Niki/MacroFold/folds/fold{}.txt".format(i+1))
        #print("pairs", pair_dict, "dangles", dangle_list)

        # Run the traceback and print the energy
        output0 = traceZ3d(0)
        #print("traceZ3d(0) is:", output0)
        output2 = -1/beta * math.log(output0)
        if output2 < 0:
            print("output2", output2)
        total_Z += output0
        energies.append(output0)
        a = round(output2, 1)
        #print("Macro:", a)

        # Get the energy from UNAfold
        output1 = subprocess.getoutput("/Users/niki/Google\ Drive/Thesis_Niki/old/unafold-3.9/src/ct-energy --suffix=DAT /Users/niki/Google\ Drive/Thesis_Niki/MacroFold/folds/fold{}.txt".format(i+1))
        b = output1 if str(output1).strip() == "+inf" else eval(output1) 
        try:
            b = float(b)
        except ValueError:
            pass  # it was a string, not a float

        #print("UNA:  ", b)
        if type(a) != type(b) or a != b:
            differences.append(("fold{}".format(i+1), all_folds[i], "Macro: {}".format(a), "UNA: {}".format(b)))

    #print()
    
    if not zero_energies:
        if not not differences:
            print("differences:")
            for e in differences:
                print(e)
    
    #for i in range(len(all_folds)):
        #print(olds)
        #print(all_folds[i], round(-1/beta*math.log(energies[i]), 1)+0, "->", energies[i])
    #print("     01234567890123456789")
    #print("Seq: {}".format(olds))
    #print("total:", -1/beta*math.log(total_G))
    #print("G =    {}".format(-1/beta*math.log(Z3[0])))
    diff = abs(total_Z - Z3[0])
    print("total_num", total_Z)
    print("we count", Z3[0])
    if abs(diff) < 1e-8:
        print("OK!")
    else:
        print(" DIFF!", diff)
        print(" SEQ", olds)
        print("total_num", total_Z)
        print("we count ", Z3[0])
        raise Exception("different")

def lock_n_load(sequence, zero_energies=False, wrap=True, no_dangles=False, likely_pairs=False, force=False):
    global Z3
    global Z5
    global N
    global s
    global olds
    global Z1
    global Z2
    global Zb

    # Call the MacroFold function and print out energy of ends pairing
    loadSequence(sequence)
    fillEnergyTables(no_dangles=no_dangles, zero_energies=zero_energies, force=force)
    fillAllowedPairsList(likely_pairs=likely_pairs)
    filleIndex()
    fillMatrices1()
    fillZ3()
    if wrap:
        fillZ5()
        fillMatrices2()

# Find differences in the number of structures found by two methods:
# 1) enumerate (exponentialy) all the different structures
# 2) set all the energies to 0 and then calculate the 
#    the partition function
# Note that this requires all the energies to be changed to 0
# (both those found in the tables AND those in this file)
def find_mismatches(num, length):
    for i in range(num):
        seq = generate_sequence(length)

        print(i, seq)
        lock_n_load(seq)

        num_configs = get_structs(seq, write=False)[0]
        configs = get_structs(seq, write=False)[2]

        all_configs = []
        counter = 0
        while len(all_configs) < num_configs:
            counter += 1
            if counter > 0 and counter % 1000 == 0:
                print("tried", counter)
            if counter > 100000:
                break
                
            lpairs, dangles = rand_sample(1)
            dangles = sorted(dangles[0])
            pairs = dict(lpairs[0])
            for base in lpairs[0]:
                pairs[pairs[base]] = base
            if not (pairs, dangles) in all_configs:
                all_configs.append((pairs, dangles))

        bad = False
        for config in all_configs:
            if config not in configs:
                bad = True
                print(config, "in Z but not listed")
        for config in configs:
            if config not in all_configs:
                bad = True
                print(config, "listed but not in Z")
        if bad:
            raise Exception("mismatch")

# Take a random sample
def rand_sample(num=1):
    global pair_dict
    global dangle_list
    sample_pairs = []
    sample_dangles = []
    for i in range(num):
        pair_dict = {}
        dangle_list = []
        traceZ3(0)
        sample_pairs.append(pair_dict)
        sample_dangles.append(dangle_list)
    return (sample_pairs, sample_dangles)

# Get the probs from the wraparound and from the actual enumeration, and compare them
def compare_probs(seq):

    lock_n_load(seq, zero_energies=True)

    num_configs, folds, configs = get_structs(seq, write=False)
    #for q,fold in enumerate(folds): print(q,fold)
    for i in range(len(seq)-4):
        for j in range(i+4, len(seq)):
            #print("looking at {},{}".format(i, j))
            calculated = (Zb[i,j]*Zb[j,i])/Z3[0]
            count = 0
            flist = []
            for d,config in enumerate(configs):
                config = config[0]
                if (i, j) in config.items():
                    count += 1
                    flist.append(folds[d])
            count = count/num_configs
            if abs(calculated-count) > 1e-11:
                print(seq)
                print(' '*i+'^'+' '*(j-i-1)+'^')
                print("i", i, "j", j)
                print(Zb[i,j], Zb[j,i])
                print("calculated", calculated*num_configs)
                print("counted   ", count*num_configs)
                print(seq)
                for e in flist:
                    print(e)
                print()
                raise Exception

# Compare probabilities.
def compare_many_probs(how_many=500, length=17):
    for i in range(how_many):
        seq = generate_sequence(length)
        if i % 10 == 0:
            print(i, seq)
        compare_probs(seq)

# Return the number of possible folds
# a given sequence can assume
def count_structures(seq, no_dangles=False, force=True):
    seq = seq.upper()
    lock_n_load(seq, zero_energies=True, wrap=False, no_dangles=no_dangles, force=force)
    #print(Z3)
    return(int(Z3[0]))

# This is for getting a set of x and y points
# corresponding to:
# x: length of randomly generated sequence
# y: total number of configurations
#import matplotlib.pyplot as plt
def save_structure_count_list(number=0, no_dangles=False):
    x = []
    y = []
    for i in range(20, 300, 40):
        with open("struct_progress.txt", 'w') as afile:
            afile.write(str(i)+'\n')
        for j in range(10):
            x.append(i)
            y.append(count_structures(generate_sequence(i), no_dangles))
    with open("x{}.txt".format(number), 'w') as myfile:
        for point in x:
            myfile.write(str(point)+'\n')
    with open("y{}.txt".format(number), 'w') as myfile:
        for point in y:
            myfile.write(str(point)+'\n')

# This is for generating the thermal ensemble
# of a given sequence and saving it as a
# .txt of points
def get_thermal_ensemble(filename="myhist.txt", seq=generate_sequence(14), num_samples=5000, 
        no_dangles=False, zero_energies=False, struct_only=False, write=False):
    # Save the sequence in the file
    if write:
        with open("./folds/structure1.txt", 'w') as mfile:
            mfile.write(str(seq)+'\n')
    #uniform = True
    global dangle_list
    global pair_dict
    energies = []
    #filename = "len_30_weighted.txt"
    #seq = 'AGCCUAUAACGAUAGGGAUUAUCCUAGCAG'
    if not zero_energies:
        lock_n_load(seq, zero_energies=zero_energies, wrap=False, no_dangles=no_dangles)
    for num, i in enumerate(range(num_samples)):
        if num % 100 == 0:
            print(num)
        dangle_list = []
        pair_dict = {}
        if zero_energies:
            lock_n_load(seq, zero_energies=zero_energies, wrap=False, no_dangles=no_dangles)
        traceZ3(0)
        if write:
            with open("./folds/structure1.txt", 'a') as mfile:
                mfile.write(str(save_structure())+'\n')
        if struct_only:
            continue
        #print(pair_dict)
        #print(dangle_list)
        # Reset the zero so we can calculate
        # the actual energy
        if zero_energies:
            lock_n_load(seq, zero_energies=False)
        energy = traceZ3d(0)
        energies.append((-1/beta)*np.log(energy))
        #print((-1/beta)*np.log(energy))

    with open(filename, 'w') as mfile:
        mfile.write(str(seq))
        mfile.write('\n')
        mfile.write(str(energies))
    #pyplot.hist(energies)
    #pyplot.show()

# How many folds are there for a given sequence length, and what 
# is the distribution of those folds?
# This function answers that question by generating
# num_samples different sequnces of length "length"
# and saves the number of possible folds for each of these
# sequences into filename, along with the length at 
# the top of the file
def folds_given_length(filename='myfolds.txt', length=30, num_samples = 10000, track_progress=False):
    numbers = []
    for i in range(num_samples):
        if track_progress and (i % 50 == 0):
            print(i)
            with open("progress.txt", 'w') as afile:
                afile.write(str(i)+'\n')
        seq = generate_sequence(length)
        num = count_structures(seq)
        numbers.append(num)

    with open(filename, 'w') as mfile:
        mfile.write("length: {}".format(length)+'\n')
        mfile.write(str(numbers))

# Save a specified number of random samples from a given sequence.
# Choose whether or not to use uniform distribution and/or dangles.
# Saves .ct files in a user-specified directory
def get_rand_samples(num_samples=1, seq=generate_sequence(30),
        zero_energies=False, no_dangles=False, output_dir="./folds/", name="sample"):
    global dangle_list
    global pair_dict
    for i in range(num_samples):
        if num % 100 == 0:
            print(num)
        dangle_list = []
        pair_dict = {}
        lock_n_load(seq, zero_energies=zero_energies, wrap=False, no_dangles=no_dangles)
        # Take the sample
        traceZ3(0) # Pairs and dangles are now stored in the pair_dict and dangle_list

        # Make a dot-bracket version of the fold
        # We'll then pass it to the make_ct function
        # to save the .ct file.
        # (do we want to put the energy anywhere?)

def recalculate_Z(threshold=1e-8, wrap=False):
    likely_pairs = [[] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if Zb[i,j]*Zb[j,i]/Z3[0] > threshold:
                likely_pairs[i].append(j)
    #print(likely_pairs)
    lock_n_load(olds, likely_pairs=likely_pairs, wrap=wrap) # only need wrap for calculating probabilites

if __name__ == "__main__":

    "testing procedure below"
    #run_tests()
    #for i in range(1000):
        #print(i)
        #seq = generate_sequence(18)
        #get_microstates(seq)
    #compare_many_probs()

    "production code to follow"
    print(count_structures(sys.argv[1], no_dangles=True), "microstates without dangles")
    print(count_structures(sys.argv[1], no_dangles=False), "microstates with dangles")

