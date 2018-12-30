# (c) 2017 Nikolaus Howe, Williams College
#
# A program which implements the recursion
# relations developed by Hofacker and Turner
# in their combinatorics papers

#import sys

global struct_dict
struct_dict = {}

def main():
    #myset = [i+1 for i in range(1000)]
    myset = [i+1 for i in range(30)]
    for i in myset:
        print(i)
    #print()
    for i in myset:
    #for i in range(30):
        #print(i)
        print(num_structs(i,1, dang=True))

def num_structs(n, eta=1, dang=False):
    global struct_dict
    if n in struct_dict:
        total = struct_dict[n]
    elif n <= 4:
        total = 1.0
    else:
        total = 0
        # 1 - no dangles
        for k in range(n-4):
            total += num_structs(k)*num_structs(n-k-2)
        # 2-4
        '''
        if dang:
            #for k in range(n-4):
                #total += num_structs(k)*num_structs(n-k-3)
            #for k in range(n-4):
                #total += num_structs(k)*num_structs(n-k-3)
            #for k in range(n-4):
                #total += num_structs(k)*num_structs(n-k-4)*(1-eta)
            # 5-8
            for k in range(1,n-4):
                total += num_structs(k-1)*num_structs(n-k-2)
            #for k in range(1,n-4):
                #total += num_structs(k-1)*num_structs(n-k-3)
            #for k in range(1,n-4):
                #total += num_structs(k-1)*num_structs(n-k-3)
            #for k in range(1,n-4):
                #total += num_structs(k-1)*num_structs(n-k-4)*(1-eta)
            # 9-12
            for k in range(n-5):
                total += num_structs(k)*num_structs(n-k-3)
            #for k in range(n-5):
                #total += num_structs(k)*num_structs(n-k-4)
            #for k in range(n-5):
                #total += num_structs(k)*num_structs(n-k-4)
            #for k in range(n-5):
                #total += num_structs(k)*num_structs(n-k-5)*(1-eta)
            # 13-16
            for k in range(1,n-5):
                total += num_structs(k-1)*num_structs(n-k-3)*(1-eta)
            #for k in range(1,n-5):
                #total += num_structs(k-1)*num_structs(n-k-4)*(1-eta)
            #for k in range(1,n-5):
                #total += num_structs(k-1)*num_structs(n-k-4)*(1-eta)
            #for k in range(1,n-5):
                #total += num_structs(k-1)*num_structs(n-k-5)*(1-eta)
        '''
        

        ''' WHERE ARE WE SUPPOSED TO MULTIPLY BY 3/8?
            LOOK IN THE PAPER '''
        # I think it's here
        total *= eta
        total += num_structs(n-1) 
        
    if not n in struct_dict:
        struct_dict[n] = total
    return total

main()
