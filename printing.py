# Nikolaus Howe, July 2016

##########################################################
# A simple program to print a nice array to the terminal #
##########################################################

from prettytable import PrettyTable

def printList(matrix, r, N):
    Elist = matrix.tolist()
    for i, row in enumerate(Elist):
        Elist[i] = row + i
    print(Elist)
    toadd = [[i for i in range(len(r))] + ['*']]
    Elist = toadd + Elist
    for i in range(len(Elist)):
        for j in range(i):
            Elist[i][j] = ' '
    for i in range(1, N):
        for j in range(i, i+3):
            if j < N:
                Elist[i][j] = '\\' 
    for i in range(N):
        Elist[i+1][i] = r[i]
    for i in range(N+1):
        for j in range(N+1):
            Elist[i][j] = str(Elist[i][j]).rjust(2)
    for i in range(1, N-3):
        for j in range(i+3, N):
            if int(Elist[i][j]) > 99:
                Elist[i][j] = "inf"
    
    p = PrettyTable()
    for row in Elist:
        p.add_row(row)
    print(p.get_string(header=False, border=True))
