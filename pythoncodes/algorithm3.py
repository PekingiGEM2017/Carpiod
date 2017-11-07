global series
series = []
global part
part = [[],[1],[2],[-2,1]]
global partrev
partrev = [[],[-1],[-2],[-1,2]]
# original parts order
global O
O = []
# recombinanse order in the gene seq
global R
R = []
sss = []
# the series of gene
global G
G = []
# 1:promoter 2:terminator -:the part is reversed

# Function:form all the possible series
#          it ensures the first part is a promoter
# Input:curN-current number
#       totN-the number of recombinanse
#       curS-current series
# Output:there is no output,but the function can change the 'series' into all the possible series
#        the lenth of series will be 2*6^(totN-1)

def formseries(curN,totN,curS):
    if curN == 1:
        formseries(curN+1,totN,curS+[1])
        formseries(curN+1,totN,curS+[-1])
    elif curN == totN+1:
        series.append(curS);
        return
    else:
        for i in range(1,4):
            formseries(curN+1,totN,curS+[i])
            formseries(curN+1,totN,curS+[-i])


# Function: find working series and not-working series in series
# Input: S-a list contains all the series
# Output:s[[s0],[s1]]
#        s0-not work
#        s1-work
def findStype(S):
    s = [[],[]]
    
    for i in S:
        temps = []
        for j in i:
            if j > 0:
                temps = temps+part[j]
            else:
                temps = temps+partrev[abs(j)]
        Ltemps = len(temps)
        flag = 0
        for j in range(Ltemps-1,-1,-1):
            if temps[j] == 1:
                flag = 1
                s[1].append(i)
                break
            if temps[j] == 2:
                flag = 1
                s[0].append(i)
                break
        if flag == 0:
            s[0].append(i)

    return s

# Function: find series can be used in the design
#           it requeires that the series can be find in another s series after reversing
# Input:S[[s0],[s1]]-which is the output of findStype(S)
#       n-the number of recombinanse
# Output:avas[[avas0],[avas1]]-available series in the design
def findtrans(S,n):
     s0 = S[0]
     s1 = S[1]
     avas = [[],[]]
     for i in s0:
         tempava = []
         for j in range(n):
             i[j] = -i[j]
             if i in s1:
                 i[j] = -i[j]
                 avas[0].append(i)
                 break
             i[j] = -i[j]
     for i in s1:
         tempava = []
         for j in range(n):
             i[j] = -i[j]
             if i in s0:
                 i[j] = -i[j]
                 avas[1].append(i)
                 break
             i[j] = -i[j]       
     return avas

# Function: find the index
# Input:S-available series
#       n-the number of recombinanse
# Output:ind[[ind0],[ind1]]-the index to find the series number ofter a part of the series is reversed
def findindex(S,n):
    s0 = S[0]
    s1 = S[1]
    ind = [[],[]]
    for i in s0:
        tempind = []
        for j in range(n):
            i[j] = -i[j]
            if i in s1:
                tempind.append(s1.index(i))
            else:
                tempind.append(-1)
            i[j] = -i[j]
        ind[0].append(tempind)
    for i in s1:
        tempind = []
        for j in range(n):
            i[j] = -i[j]
            if i in s0:
                tempind.append(s0.index(i))
            else:
                tempind.append(-1)
            i[j] = -i[j]
        ind[1].append(tempind)             
    return ind

# Function: inorder t
# Input: statenum = len(rec)+1 = len(res)
#s-[1,statenum]

def way2_1v2F(rec,res,s,currentS,originalS,RecS,statenum):
    if s == statenum:
        O.append(originalS)
        R.append(RecS.copy())
        G.append(currentS)
        return
    if res[s]^res[s - 1]:
        p = res[s-1]
        q = avas[p].index(currentS)
        tempind = ind[p][q]
        for i in range(statenum-1):
            if tempind[i] != -1:
                RecS[i] = RecS[i]+rec[s-1]

                way2_1v2F(rec,res,s+1,avas[res[s]][tempind[i]],originalS,RecS,statenum)
                RecS[i] = RecS[i][0:len(RecS[i])-1]
    else:
        way2_1v2F(rec,res,s+1,currentS,originalS,RecS,statenum)

def way2_1v2(rec,res,n):
    Lrec = len(rec)
    




NN = 2
formseries(1,NN,[])
ss = findStype(series)
global avas
avas = ss
global ind
ind = findindex(ss,NN)


#for i in avas[0]:
#    recs = []
#    for j in range(NN):
#        recs.append('')
#    way2_1v2F('AB',[0,0,1],1,i,i,recs,3)

#for i in avas[1]:
#    recs = []
#    for j in range(NN):
#        recs.append('')
#    way2_1v2F('AB',[1,1,0],1,i,i,recs,3)

OO = O.copy()
RR = R.copy()
GG = G.copy()
O = []
R = []
G = []

res = [0,0,1]
rec = 'AB'
Lres = len(res)
Lrec = len(rec)
for i in avas[res[0]]:
    recs = []
    for j in range(NN):
        recs.append('')
    way2_1v2F(rec,res,1,i,i,recs,Lres)

LO = len(O)
print(res)
for i in range(LO):
    for j in range(Lrec):
        print(R[i][j],O[i][j],R[i][j][::-1],end = ' ')
    print()
