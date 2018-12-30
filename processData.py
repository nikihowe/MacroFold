# for reading rawdata/loop.txt and outputting :
# -data/hairpin.csv
# -data/bulge.csv 
# -data/internal.csv
def convertHairpinBulgeInternal():
    # open files
    loopFile = open("rawdata/loop.txt", 'r')
    hairpinFile = open("data/hairpin.csv", 'w')
    bulgeFile = open("data/bulge.csv", 'w')
    internalFile = open("data/internal.csv", 'w')

    ## write headers
    hairpinFile.write("Length,Energy\n")
    bulgeFile.write("Length,Energy\n")
    internalFile.write("Length,Energy\n")

    skips = 4 # hardcoded number of lines to skip
    for line in loopFile:
        if(skips > 0):
            skips -= 1
            continue
        tokens = line.split()
        
        # internal energy column
        if tokens[1] is not ".":
            internalFile.write(tokens[0])
            internalFile.write(",")
            internalFile.write(tokens[1])
            internalFile.write("\n")

        # bulge energy column
        if tokens[2] is not ".":
            bulgeFile.write(tokens[0])
            bulgeFile.write(",")
            bulgeFile.write(tokens[2])
            bulgeFile.write("\n")

        # hairpin energy column
        if tokens[3] is not ".":
            hairpinFile.write(tokens[0])
            hairpinFile.write(",")
            hairpinFile.write(tokens[3])
            hairpinFile.write("\n")
    
    loopFile.close()
    hairpinFile.close()
    bulgeFile.close()
    internalFile.close()

def intToBase(num):
    if num == 0: 
        return 'A'
    elif num == 1:
        return 'C'
    elif num == 2:
        return 'G'
    elif num == 3:
        return 'U'
    else:
        print("gaah!")
        return 'N'

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
    else:
        return ('N', 'N')

def intToTwoBase(num):
    first  = intToBase(num // 4)
    second = intToBase(num %  4)
    return (first, second)

def convertStack():
    stackFile = open("rawdata/stack.txt", 'r')
    outputFile = open("data/stack.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    skip = 26
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
            skip = 11
    stackFile.close()
    outputFile.close()

def convertTStackH():
    stackFile = open("UNAdata/tstackh.DAT", 'r')
    outputFile = open("data/tstackh.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy == "." or energy == "inf":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
    stackFile.close()
    outputFile.close()

def convertTStackI():
    stackFile = open("UNAdata/tstacki.DAT", 'r')
    outputFile = open("data/tstacki.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy == "." or energy == "inf":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
    stackFile.close()
    outputFile.close()

def convertTStackE():
    stackFile = open("UNAdata/tstacke.DAT", 'r')
    outputFile = open("data/tstacke.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy == "." or energy == "inf":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
    stackFile.close()
    outputFile.close()

def convertTStackM():
    stackFile = open("UNAdata/tstackm.DAT", 'r')
    outputFile = open("data/tstackm.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy == "." or energy == "inf":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
    stackFile.close()
    outputFile.close()

def convertDangles():
    dangleFile = open("rawdata/dangle.txt", 'r')
    dangle5 = open("data/dangle5.csv", 'w')
    dangle3 = open("data/dangle3.csv", 'w')

    dangle5.write("Dangle, i, j, Energy\n")
    dangle3.write("i, j, Dangle, Energy\n")
    skip = 10
    threePrime = True
    index5prime = 0

    for line in dangleFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue
                #energy = "Inf"
            dangleBase = intToBase(i % 4)
            base3prime = intToBase(i // 4)
            base5prime = intToBase(index5prime)
            if threePrime:
                dangle3.write("{0},{1},{2},{3}\n".format(base3prime, base5prime, dangleBase, energy))
            else:
                dangle5.write("{0},{1},{2},{3}\n".format(dangleBase, base3prime, base5prime, energy))
            
        index5prime += 1
        skip = 10
        if index5prime >= 4:
            threePrime = False
            index5prime = 0
        
    dangleFile.close()
    dangle5.close()
    dangle3.close()
    

def convertTerminalStack():
    # open files
    inputFile = open("rawdata/tstack.txt", 'r')
    outputFile = open("data/tstack.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    skip = 26
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in inputFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue;
                #energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
            skip = 11
    inputFile.close()
    outputFile.close()
    

# parse the 1x1 energy file into .csv
def convert1x1():
    # open the files to read and write from/to
    file1x1 = open("rawdata/int11.txt", 'r')
    outputFile = open("data/int11.csv", 'w')

    outputFile.write("5primeOuter,5primeMismatch,5primeInner,3primeInner,3primeMismatch,3primeOuter,Energy\n")

    skip = 30
    mismatch5primeIndex = 0
    outerIndex = 0
    for line in file1x1:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue
                ##energy = "Inf"
            mismatch3prime = intToBase(i % 4)
            (inner5prime, inner3prime) = intToPair(i // 4)
            mismatch5prime = intToBase(mismatch5primeIndex)
            (outer5prime, outer3prime) = intToPair(outerIndex)
            outputFile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(outer5prime, mismatch5prime, inner5prime, inner3prime, mismatch3prime, outer3prime, energy))
        mismatch5primeIndex += 1
        if mismatch5primeIndex >= 4:
            mismatch5primeIndex = 0
            outerIndex += 1
            skip = 11
    file1x1.close()
    outputFile.close()


# parse the 2x1 energy file into .csv
def convert2x1():
    # open the files to read and write from/to
    file2x1 = open("rawdata/int21.txt", 'r')
    outputFile = open("data/int21.csv", 'w')

    outputFile.write("5primeOuter,5primeMismatch,5primeInner,3primeInner,3primeInnerMismatch,3primeOuterMismatch,3primeOuter,Energy\n")

    skip = 30
    mismatch5primeIndex = 0
    outerIndex = 0
    for n,line in enumerate(file2x1):
        #print(n+1)
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            #print(tokens)
            energy = tokens[i]
            #print(energy)
            #print(i)
            #print(outerIndex)
            if energy is ".":
                continue
                ##energy = "Inf"
            #print(outerIndex)
            (outer5prime, outer3prime) = intToPair(outerIndex // 4)
            (inner5prime, inner3prime) = intToPair(i // 4)
            mismatch5prime = intToBase(mismatch5primeIndex)
            outermismatch3prime = intToBase(i % 4)
            innermismatch3prime = intToBase(outerIndex % 4)
            outputFile.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(outer5prime, mismatch5prime, inner5prime, inner3prime, innermismatch3prime, outermismatch3prime, outer3prime, energy))
        mismatch5primeIndex += 1
        if mismatch5primeIndex >= 4:
            mismatch5primeIndex = 0
            outerIndex += 1
            skip = 11
    file2x1.close()
    outputFile.close()


# parse the 2x1 energy file into .csv
def convert2x2():
    # open the files to read and write from/to
    file2x2 = open("rawdata/int22.txt", 'r')
    outputFile = open("data/int22.csv", 'w')

    outputFile.write("5primeOuter,5primeOuterMismatch,5primeInnerMismatch,5primeInner,3primeInner,3primeInnerMismatch,3primeOuterMismatch,3primeOuter,Energy\n")

    skip = 41
    row = 0
    outerIndex = 0
    for n,line in enumerate(file2x2):
        #print(n+1)
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for col in range(len(tokens)):
            #print(tokens)
            energy = tokens[col]
            #print(energy)
            #print(i)
            #print(outerIndex)
            if energy is ".":
                continue
                ##energy = "Inf"
            #print(outerIndex)
            (outer5prime, outer3prime) = intToPair(outerIndex // 6)
            (inner5prime, inner3prime) = intToPair(outerIndex %  6)
            (outermismatch5prime, outermismatch3prime) = intToTwoBase(row)
            (innermismatch5prime, innermismatch3prime) = intToTwoBase(col)
            outputFile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(outer5prime, outermismatch5prime, innermismatch5prime, inner5prime, inner3prime, innermismatch3prime, outermismatch3prime, outer3prime, energy))
        row += 1
        if row >= 16:
            row = 0
            outerIndex += 1
            skip = 10
    file2x2.close()
    outputFile.close()

# Convert Example
def convInt11():
    stackFile = open("UNAdata/sint2.DAT", 'r')
    outputFile = open("data/UNAint11.csv", 'w')

    outputFile.write("5primeOuter,5primemismatch,5primeInner,3primeInner,3primemismatch,3primeOuter,Energy\n")
    outer5primeIndex = 0
    inner5primeIndex = 0

    # Get all the energies into a list
    allenergies = []
    for line in stackFile:
        for token in line.split():
            allenergies.append(token)

    print(allenergies)
    # Follow the parsing procedure that they do
    # in energy.c
    index = 0
    for b in range(6):
        for i in range(5):
            for c in range(6):
                for j in range(5):
                    # Get the energy and get ready for the next one
                    energy = allenergies[index]
                    index += 1

                    if energy == "inf":
                        continue
                        ##energy = "Inf"

                    # Get the indices
                    outer5prime, outer3prime = intToPair(b)
                    inner5prime, inner3prime = intToPair(c)
                    mismatch5prime = intToBase(i)
                    mismatch3prime = intToBase(j)
                    print("sequence: {}{}{}{}{}{}".format(outer5prime, mismatch5prime, inner5prime, inner3prime, mismatch3prime, outer3prime))
                    print("energy: %s"%energy)
                    if (b == 6 or c == 6):
                        "do nothing"
                        #outputFile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(outer5prime, 5primemismatch, inner5prime, inner3prime, 3primemismatch, outer3prime, "inf"))
                    elif (i == 4 or j == 4):
                        "do nothing"
                        #outputFile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(outer5prime, 5primemismatch, inner5prime, inner3prime, 3primemismatch, outer3prime, 0))
                    else:
                        outputFile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(outer5prime, mismatch5prime, inner5prime, inner3prime, mismatch3prime, outer3prime, energy))

    stackFile.close()
    outputFile.close()


#convertHairpinBulgeInternal()
#convertTerminalStack()
#convertDangles()
#convertStack()
#convert1x1()
#convert2x1()
#convert2x2()
#convertTStackH()
#convertTStackE()
#convertTStackI()
#convInt11()
convertTStackM()
