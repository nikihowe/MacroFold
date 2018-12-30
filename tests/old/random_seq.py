import sys
import random

length = int(sys.argv[1])

print(''.join([['A', 'C', 'G', 'U'][random.randrange(4)] for i in range(length)]))

