import math 

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
    
if __name__ == '__main__':

    numbers = [1,2,3,4]

    res = permutations(numbers)

    with open('m.txt', 'w') as file:
        file.write('\n')

        file.write(str(len(res)))
        
        file.write('\n')
        for perm in res:
            file.write(' '.join(map(str,perm)))
            file.write('\n')