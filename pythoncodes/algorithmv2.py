#Input
# rec is the order of recombinanse, which is equal to the number of state - 1
# the length of rec and seq should be the same
# res is the result of output which is a series of number(1 and 0)
# Outputs
# number 1 represents a promoter, number 2 represents a terminator
# the negative number represents the part is reversed
# [1 -2] means work; [-1 -2] and[1 2]means not work

def Way2_1(rec, res):
    print(res)
    print(rec)
    if res[0] == 0:
        Way2_1F(1,[1,2],[1,2],['',''])
        Way2_1F(1,[-1,-2],[-1,-2],['',''])
    elif res[0] == 1:
        Way2_1F(1,[1,-2],[1,-2],['',''])
    case = len(O)
    for i in range(case):
        print(R[i][0],O[i][0],R[i][0][::-1],R[i][1],O[i][1],R[i][1][::-1])

def Way2_1F(s, currentS, originalS, RecS):
    if s == stanum:
        O.append(originalS)
        R.append(RecS)
        return
    if res[s]^res[s - 1]:
        if currentS == [1,-2]:
            Way2_1F(s + 1, [1,2], originalS, [RecS[0],RecS[1]+rec[s - 1]])
            Way2_1F(s + 1, [-1,-2], originalS, [RecS[0]+rec[s - 1],RecS[1]])                    
        elif currentS == [1,2]:
            Way2_1F(s + 1, [1,-2], originalS, [RecS[0],RecS[1]+rec[s - 1]])
        elif currentS == [-1,-2]:
            Way2_1F(s + 1, [1,-2], originalS, [RecS[0]+rec[s - 1],RecS[1]])
    else:
        Way2_1F(s + 1, currentS, originalS, RecS)

# original parts order
global O
O = []
# recombinanse order in the gene seq
global R
R = []
# the recombinanse input order
global rec
#rec = input('please input the order of recombinanse')
rec = 'AB'
global res
#res = input('plese input the order of outputs')
res = [1,0,0]
global recnum
stanum = len(res)
Way2_1(rec,res)


    
    
