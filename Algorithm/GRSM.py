#  Reverse a part
def rev(oblock): 
    block = oblock.copy()
    block.reverse()
    block = [int(i)*(-1) for i in block]
    return block

# 将输入的Inputnumber中parts的序号对应成相应的元件
def initinputs(inputnumber, parts): 
    inputs = []
    for i in inputnumber:
        temp = parts[abs(i)].copy()
        if i < 0:
            temp = rev(temp)
        inputs.append(temp)
    cnt = 0
    for i in inputs:
        n = len(i)
        for j in range(n):
            if i[j] == 3:
                i[j]=i[j]+cnt
                cnt+=1
            elif i[j] == -3:
                i[j]=i[j]-cnt
                cnt+=1
    # 最后一个元素为初始线路基因的数量
    inputs.append([cnt])
    return inputs

# 根据状态得出每种状态下元件的链接方式
def change(inputs, state):
    outputs = []
    for i in state:
        output = []
        for j in i:
            j = int(j)
            temp = inputs[abs(j) - 1].copy()   
            if j < 0:
                temp = rev(temp) 
            output = output + temp
        outputs.append(output)
    return outputs

# 根据元件的连接方式得到表达的基因
def outG(outputs):
    outgenes = []
    flag = 0
    for i in outputs:
        outgene = []
        n = len(i)
        for j in range(n):
            if i[j] > 2:
                cnt = j
                while cnt > 0:
                    cnt = cnt-1
                    if i[cnt]==2 or abs(i[cnt])>2:
                        break
                    elif i[cnt] == 1:
                        num = abs(i[j]) - 2
                        outgene.append(num)
                        flag = flag + 1
                        break
            if i[j] < -2:
                cnt = j
                while cnt < n-1:
                    cnt = cnt+1
                    if i[cnt]==-2 or abs(i[cnt])>2:
                        break
                    elif i[cnt] == -1:
                        num = abs(i[j]) - 2
                        outgene.append(num)
                        flag = flag + 1                
                        break
                    # flag为表达的基因数量
        outgenes.append(outgene)
    # 保证状态0-1 0-2 1-3 2-4不同
    s1 = [0,0,1,2]
    s2 = [1,2,3,4]
    for i in range(4):
        if outgenes[s1[i]] == outgenes[s2[i]]: 
            flag = 0
            break
    outgenes.append([flag])
    return outgenes

# 转化为01矩阵
def outmat(outg):
    arr = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
    for i in range(5):
        for j in outg[i]:
            arr[i][j-1] = 1
    arr1 = []
    for i in range(5):
        temp = 0
        for j in range(5):
            if arr[j][i]==1:
                temp = temp+2**(4-j)
        arr1.append(temp)
    arr1.sort()
    arr1.reverse()
    arr2 = [[],[],[],[],[]]
    for i in range(5):
        for j in range(5):
            arr2[4-j].append(arr1[i]%2)
            arr1[i] = arr1[i]//2
    arr3 = []
    for i in arr2:
        temp = 0
        for j in range(5):
            if i[j] == 1:
                temp = temp+2**j
        arr3.append(temp)
    if arr3[1] > arr3[2]:
        arr3[1],arr3[2] = arr3[2],arr3[1]
        arr3[3],arr3[4] = arr3[4],arr3[3]
    return arr3

# 元件初始化
# 注意此部分在每个part里面能够受到控制的基因最多为两个，两端的part受控制的基因为1个
def initparts():
    P = 1  # promoter
    T = 2  # Terminator
    G = 3  # Gene with bi-directional terminator
    
    parts = ([],)
    parts = parts+([G],)
    parts = parts+([G,P],)
    parts = parts+([T,T*(-1),P],)
    parts = parts+([T*(-1),T],)
    parts = parts+([],)
    parts = parts+([P*(-1),P],)
    parts = parts+([T],)
    parts = parts+([G,G*(-1)],)
    parts = parts+([G*(-1),P],)
    parts = parts+([P],)
    parts = parts+([P*(-1),G,P],)
    parts = parts+([G,G*(-1),P],)
    parts = parts+([P*(-1),G,G*(-1),P],)
    parts = parts+([T*(-1),P],)
    parts = parts+([T,T*(-1)],)
    parts = parts+([G,P,G],)
    parts = parts+([G,P,G,P],)
    parts = parts+([P,G,P],)
    parts = parts+([P,G],)
    parts = parts+([P*(-1),P,G,P],)
    parts = parts+([G,P,G,G*(-1)],)
    parts = parts+([P,G,G*(-1),P],)
    parts = parts+([P*(-1),G,P,G,P],)
    parts = parts+([G,P,G,G*(-1),P],)
    parts = parts+([P*(-1),G,P,G,G*(-1),P],)
    return parts

# 递归输出文件
# 格式为：
# 第一行 DNA链的parts编码序列
# 第二行 在每种状态下表达的基因，从编号3开始递增，在每一个输出组合里面不会有重复的基因编号
# 第三行 空行
# 注意：修改保存路径或者文件名，以及换行符在不同系统下可能不同
def outdata(Inputnumber, part, statenumber, n, partnumber):
    if n == part:             
        global Parts
        i = initinputs(Inputnumber,Parts)
        cnt = i[len(Inputnumber)][0]
        # cnt为初始序列携带基因的数量，如果多于5则舍弃
        if cnt > 5:
            return
        i.pop()
        global State
        o = change(i,State)
        g = outG(o)
        # 保证表达基因的数量和起始序列基因的数量一致
        if(g[statenumber][0] == cnt):
            g.pop()
            fo = outmat(g)
            sI = str(Inputnumber)
            sg = str(fo)
            f = open('text.txt','a')
            f.write('\n'+ '\n'+sI)
            f.write('\n'+ sg)
            f.close()
        return
    for i in range((-1)*partnumber,partnumber+1):
        # 0和空不需要翻转
        if i == 0 or i == -5:
            continue
        Inputnumber[n] = i
        outdata(Inputnumber, part, statenumber, n+1,partnumber)
        
# 不同识别序列下的状态，记得修改State，共有16种0.0    
global State
State = ([1,2,3,4,5,6,7],[1,2,3,6,7],[1,-6,-5,3,4,-2,7],[1,-6,-3,-2,7],[1,-6,-3,5,4,-2,7])

L = len(State)  # 状态数量

Inputnumber = [0,0,0,0,0,0,0]

global Parts
Parts = initparts()

outdata(Inputnumber,7, L,0,25)




