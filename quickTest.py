#!/usr/bin/env python
from  findParent import countSplit
from  findParent import countDup
from  findParent import reductionSubset
from  findParent import reductionCount
from  findParent import findSetInitial_GG
from  findParent import findSetInitial_SG
from  findParent import findSetInitial_SS

# quick test:
XC = 'hj|fg|e'
XA = 'fg|hj|e'
SO = 'fj|fg|h|i|e'
PA = 'fhj|eg|k|i'
PS ='fj|h|i|g|e'
PP='fghijkabcde'
VP='fj|k|h|g|e'
VV='fj|hg|e'
YP='fjh|g|e'
SE='fj|k|hg|e'
SF='fjhg|e'
EC='abcdefghijk'
A = findSetInitial_GG(XC,XA,2,2)
B = findSetInitial_GG(PP,PS,0,4)
C = findSetInitial_GG(VP,VV,4,2)
D = findSetInitial_GG(EC,SF,0,1)

E = findSetInitial_SG(A,SO,4)
F = findSetInitial_SG(B,PA,3)
G = findSetInitial_SG(D,SE,3)

H = findSetInitial_SS(E,F)
I = findSetInitial_SG(G, YP,2)

J= findSetInitial_SS(C,I)
K= findSetInitial_SS(J,H)

L= findSetInitial_SS(J,K)
print('A:',A)
print('B:',B)
print('C:',C)
print('D:',D)
print('E:',E)
print('F:',F)
print('G:',G)
print('H:',H)
print('I:',I)
print('J:',J)
print('K:',K)
print('L:',L)
