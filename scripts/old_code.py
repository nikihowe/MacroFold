
# The below method adds in every combination of stacked dangles
# adjacent to every normal paired stack. Then I attempt to remove
# all the cases which are not allowed. This proved to be really
# difficult. The new approach now adds in the dangles only
# in the two places they are allowed to appear:
# in multibranches and the external loop.

"""
# Now we go through the list and, for every space next to every
# paren, we try replacing it with a < or a >, where 
# the stacking occurs between the < and the thing 
# to its right, or between the > and the thing 
# to its left
old_folds = []
while old_folds != all_folds:
    old_folds = all_folds
    for fold in all_folds:
        for i, char in enumerate(fold):
            if char == '.':
                if i > 0 and (fold[i-1] == '(' or fold[i-1] == ')'):
                    # now we need to append the fold in which this 
                    # base is now stacked to the left. We'll do a 
                    # similar one for stacking on the right.
                    # We should repeat the whole procedure on all_folds
                    # until there is no further updating
                    new_fold = fold[:i] + '>' + fold[i+1:]
                    if not new_fold in all_folds:
                        all_folds.append(new_fold)
                if i < len(fold)-1 and (fold[i+1] == '(' or fold[i+1] == ')'):
                    new_fold = fold[:i] + '<' + fold[i+1:]
                    if not new_fold in all_folds:
                        all_folds.append(new_fold)

# Now we have all possible folds. Problematically, we also have 
# dangles which are not accounted for in our model,
# namely dangles inside hairpins, loops and bulges.
# We now remove these.

# Remove the dangles inside hairpins
new_folds = set()
for fold in all_folds:
    new_fold = fold
    start = -1
    for i, char in enumerate(fold):
        if char == '(':
            start = i
        if char == ')':
            if start >= 0:
                # inclusive
                hairpin = (start, i)
                new_fold = new_fold[:start+1] + '.'*(i-start-1) + new_fold[i:]
                new_folds.add(new_fold)
            start = -1
new_folds.add('.'*len(all_folds[0]))
all_folds = new_folds
print(new_folds)

# Remove the dangles inside loops and bulges
new_folds = set()
for fold in all_folds:
    new_fold = fold
    # Stage 0: start
    # Stage 1: have found external opening bracket
    # Stage 2: have found internal opening bracket
    # Stage 3: have found internal closing bracket
    # Stage 4: have found external closing bracket
    stage = 0
    for char in fold:
        if char = '(':
            if stage = 0:
                stage = 1
            if stage >= 1:
                stage 
    
    

# Finally, we need to remove tstacks where the bases
# which are stacking could also form a pair.
"""
