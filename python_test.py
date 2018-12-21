#!/usr/bin/env python3
print("Forget your Java things for now :) And use Python!, Java Java Java")

sequence = 'TCCGATTGACTTAGA'
print(sequence.index('TCTCAACT'))
print(sequence.count('TCTCAACT'))


#def find_all(a_str, sub):
#    start = 0
#    while True:
#        start = a_str.find(sub, start)
#        if start == -1: return
#        yield start
#        start += len(sub) # use start += 1 to find overlapping matches

#list(find_all('spam spam spam spam', 'spam')) # [0, 5, 10, 15]


#def find_all(a_str, sub):
#    indexes = []
#    k = 0
#    while k < len(a_str):
#        k = a_str.find(sub, k)
#        if k == -1:
#            return indexes
#        else:
#            indexes.append(k)
#            k += 1
#    return indexes


def find_all(a_str, sub):
    indexes = []
    k = 0
    while k < len(a_str):
        if a_str[k:k+len(sub)] == sub:
            indexes.append(k)
        k += 1
    return indexes
    #print(indexes)

list(find_all('ACUUGACUCUAT', 'CU'))