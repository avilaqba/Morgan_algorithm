import pandas as pd
import numpy as np
lines = open(r'morgan.mol').readlines()
Nlist = [[],[],[]]
Dict = dict()
c = 1
Hlist = []
for line in lines:
    if len(line.split()) == 7:
        Nlist[0].append(int(line.split()[0]))
        Nlist[1].append(int(line.split()[1]))
        Nlist[2].append(int(line.split()[2]))
        Hlist.append([int(line.split()[0]),int(line.split()[1])])
    if len(line.split()) == 16:
        Dict[c] = line.split()[3]
        c += 1

df1 = pd.DataFrame(
        {'Atom1' : Nlist[0],
         'Atom2' :  Nlist[1],
         'Bonds' : Nlist[2],
         }
     )
df2 = pd.DataFrame(
        {'Atom1' : Nlist[1],
         'Atom2' :  Nlist[0],
         'Bonds' : Nlist[2],
         }
     )
     
df = pd.concat([df1, df2])
I = df.groupby(['Atom1']).sum()

Elist, Dlist, Ilist = [], [], []
for i in I['Bonds'].index:
    Elist.append(Dict[I['Bonds'].index[i-np.min([I.index])]])
    Dlist.append(I['Bonds'][I['Bonds'].index[i-np.min([I.index])]])
    Ilist.append(i)

Emlist = [Dlist]  
No = len(set(Emlist[0]))
em = 1         
NoUp = 1000000
while NoUp > No:
    Emlist.append([])
    for e in Ilist:
        Sum = 0
        for f in Ilist:
            if ([e,f] in Hlist) or ([f,e] in Hlist):
                Sum += Emlist[em-1][f-1]
        Emlist[em].append(Sum)
    No = len(set(Emlist[em-1]))
    NoUp = len(set(Emlist[em]))
    em += 1
    
Llist = Emlist[-2]
Rlist = [Llist.index(max(Llist))+1] 
StPoint = Llist.index(max(Llist))
Mlist = []
for e in set(Ilist):
    Flist = []
    for f in set(Ilist):
        if ([e,f] in Hlist) or ([f,e] in Hlist):
            Flist.append(Llist[f-1])
    Mlist.append(Flist)

c = 0
while c < len(Ilist):
    while len(Mlist[StPoint]) != 0:
        for ee in [i for i, e in enumerate(Emlist[-2]) if e == max(Mlist[StPoint])]:
            if ee+1 not in Rlist:
                Rlist.append(ee+1)
        Mlist[StPoint].remove(max(Mlist[StPoint]))
    StPoint = Rlist[c]-1
    c += 1
print 'The script prints out order of elements according to Morgan algorithm. In order to differentiate among same elements of the compound, the row numbers of elements in the .mol file has been considered.\n'    
for r in range(len(Rlist)):    
    print 'Element %s in the %dth row' %(Dict[Rlist[r]],Rlist[r])
