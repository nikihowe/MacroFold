# Nikolaus Howe July 2016
# Make a set of pairs and save in pairs.txt

import sys

def intToBase(num):
    if num == 0: 
        return 'A'
    elif num == 1:
        return 'C'
    elif num == 2:
        return 'G'
    elif num == 3:
        return 'U'

def intToPair(num):
    if num == 0: 
        return ('A', 'U')
    elif num == 1:
        return ('C', 'G')
    elif num == 2:
        return ('G', 'C')
    elif num == 3:
        return ('U', 'A')
    elif num == 4:
        return ('G', 'U')
    elif num == 5:
        return ('U', 'G')

def intToTwoBase(num):
    first  = intToBase(num // 4)
    second = intToBase(num %  4)
    return (first, second)

if __name__=="__main__":
    output_file = open('pairs.txt', 'w')
    output_list = []
    for num in range(16):
        output_list.append(intToTwoBase(num))
    for pair in output_list:
          output_file.write("%s\n" % str(pair))
