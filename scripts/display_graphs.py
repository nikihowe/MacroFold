# A quick script to display plots
# Nikoalus Howe, Jan 9 2017

from matplotlib import pyplot
import numpy as np
from scipy import stats

with open("file.txt") as mfile:
    a = mfile.read()

a = a.strip()
a = a.split('\n')

x = eval(a[0])
y = eval(a[1])
print(len(y))

for i in range(0, len(y), 2000):
    data = y[i:i+2000]
    print(i, i+2000)
    print(" mean:  ", np.mean(data))
    print(" median:", np.median(data))
    print(" mode:  ", stats.mode(data)[0][0])
    pyplot.hist(data, bins=25)
    pyplot.title("len = {}".format(x[i]))
    pyplot.savefig('len{}size2000bins25.png'.format(i))
    pyplot.clf()
