# Niki Howe, Jan 2017

# This program takes a sequence and outputs 
# to the current directory
# a set of .ct files
# named fold1.ct through foldN.ct
# corresponding to all possible 
# secondary structures, dangles included.


# imports
import sys
import subprocess
import random

# Main program
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

# Get the sequence
#:seq = "CCCAAAGGG"
get_structs("CCCAAAGGG")

# Write all possible folds to the folds directory
#main(str(seq))

# Below used for generating the plots
#'''
#sequences = []
#for n in range(5, 21, 3):
    #for i in range(2000):
        #if i % 100 == 0:
            #print(n, i)
        #seq = generate_sequence(n)
        #combos = main(str(seq), True)
        #sequences.append((n, combos))
#for seq in sequences:
    #print(seq)
#seq_length = []
#num_folds = []
#for x, y in sequences:
    #seq_length.append(x)
    #num_folds.append(y)
#with open("len.txt", 'w') as out_file:
    #out_file.write(str(seq_length))
#with open("num.txt", 'w') as out_file:
    #out_file.write(str(num_folds))
#'''
