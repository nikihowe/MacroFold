# a quick test to see which way is faster in filling arrays.
# (c) Nikolaus Howe 2017

import numpy as np
import time

def calculate_filling_times(SIZE):
    output = "{} ".format(SIZE)
    a = np.array([1.0 for i in range(SIZE*SIZE)])
    a.shape = (SIZE, SIZE)

    start_time = time.time()
    for d in range(1, SIZE):
        for i in range(SIZE-d):
            j = i+d 
            a[i, j] = 0.1*a[i+1, j] + 0.1*a[i, j-1]
    end_time = time.time()
    np.set_printoptions(precision=2, suppress=True)
    output += str(end_time - start_time)
    output += ' '

    start_time = time.time()
    for j in range(1, SIZE):
        for i in reversed(range(0, j)):
            a[i, j] = 0.1*a[i+1, j] + 0.1*a[i, j-1]
    end_time = time.time()
    output += str(end_time - start_time)
    output += '\n'

    return output

if __name__ == "__main__":
    toprint = ""
    for size in [250, 500, 750, 1000, 1250, 1500, 1750]:
        for i in range(10):
            toprint += calculate_filling_times(size)
    with open("times.txt", 'w') as afile:
        afile.write(toprint)
