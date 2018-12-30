with open("tstack.csv", 'r') as mfile:
    a = mfile.read()
    a = a.strip().split('\n')

for i, line in enumerate(a):
    if i == 0:
        continue
    else:
        line = line.split(',')
        line = line[:-1]
        line.append('0')
        a[i] = line
#U,U,U,G,U,U,G,2.7
with open("tstack.csv", 'w') as mfile:
    for i,e in enumerate(a):
        print(i,e)
        for i in range(4):
            mfile.write(e[i])
            mfile.write(',')
        mfile.write(e[4])
        mfile.write('\n')
