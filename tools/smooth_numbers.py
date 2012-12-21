#!/usr/bin/python

import itertools

def smooth_numbers(N):
    fact = [2,3,5,7]
    big_helper = [1, 11, 13]
    numbers = []
    limits = [int(N**(1./x)) for x in fact]
    exponents = []
    exponents = [range(i) for i in limits]
    
    for (a,b,c,d) in [(a,b,c,d) for a in exponents[0] for b in exponents[1] for c in exponents[2] for d in exponents[3]]:
        number = (2**a)*(3**b)*(5**c)*(7**d)
        if number > N:
            break
        for i in big_helper:
            if i*number <= N:
                numbers.append(i*number)
    return numbers            

sn = smooth_numbers(2048)
print sn


