import math
from helpers import read_file

def mortal_fib(n, k):
    age = [0] * (k)
    age[0] = 1
    for day in range(1,n):
        newborns = sum(age[1:])
        for group in range(k-1, 0, -1):
            age[group] = age[group-1]
        age[0] = newborns
        
    return sum(age)
            
def fib(n, k):
    previous1, previous2 = 1, 1

    for i in range(2, n):
        current = previous1 + k * previous2
        previous2 = previous1
        previous1 = current
    return current

test = 0

def permutations(numbers):
    res = []

    if len(numbers) == 0:
        return [numbers]
    
    for i in range(0, len(numbers)):
        for perms in permutations(numbers[:i] + numbers [i+1:]):
            new = [numbers[i]] + perms
            res.append(new)
        
    return res

def subsets(elements):
    res = [[]]
    for elem in elements:
        new = []
        for elem_old in res:
            new.append(elem_old + [elem])
        res.extend(new)

    return res 

def set_operations(n, set1, set2):
    union = set1 | set2

    ns = set([i+1 for i in range(n)])

    intersection = set1 & set2


    difference12 = set1 - set2
    difference21 = set2 - set1

    set1_complement = set1 ^ ns
    set2_complement = ns ^ set2 

    return [
        union,
        intersection,
        difference12,
        difference21,
        set1_complement, 
        set2_complement 
    ]

def main():
    n = 10
    set1 = {1, 2, 3, 4, 5}
    set2 = {2, 8, 5, 10}

    print(set_operations(n,set1,set2))

if __name__ == '__main__':
    main() 
