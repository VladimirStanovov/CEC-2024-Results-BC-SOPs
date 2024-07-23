import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import math
from ranking import Ranking, FRACTIONAL

def MannWhitneyU_myZ(Sample1,Sample2):        
    NewSample = np.concatenate((Sample1,Sample2),axis=0)
    NewRanks, Groups = get_fract_ranks_and_groups(NewSample)
    SumRanks = 0
    SumRanks2 = 0
    for i in range(Sample1.shape[0]):
        SumRanks += NewRanks[i]
        SumRanks2 += NewRanks[Sample1.shape[0]+i]
    U1 = SumRanks - Sample1.shape[0]*(Sample1.shape[0]+1.0)/2.0
    U2 = SumRanks2 - Sample2.shape[0]*(Sample2.shape[0]+1.0)/2.0
    Umean = Sample1.shape[0]*Sample2.shape[0]/2.0
    GroupsSum = 0
    for index in Groups:
        GroupsSum += (index*index*index - index)/12
    N = Sample1.shape[0]+Sample2.shape[0]
    part1 = Sample1.shape[0]*Sample2.shape[0]/(N*(N-1.0))
    part2 = (N*N*N-N)/12.0
    Ucorr2 = math.sqrt(part1*(part2-GroupsSum))
    if(Ucorr2 != 0):
        Z1 = (U1 - Umean)/Ucorr2
        Z2 = (U2 - Umean)/Ucorr2
    else:
        return (0,0)
    if(Z1 <= Z2):
        if(Z1 < -2.58):
            return (-1, Z1)
    else:
        if(Z2 < -2.58):   
            return (1, Z1)
    return (0, Z1)
def ranks(ranking):
    return list(ranking.ranks())   
def get_fract_ranks(data):
    sort_index = np.argsort(-data)
    sort_list = -np.sort(-data)
    new_ranks = ranks(Ranking(sort_list, FRACTIONAL))
    index_rank = np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        new_rank_inv = data.shape[0] - new_ranks[i] - 1
        index_rank[sort_index[i]] = new_rank_inv
    return index_rank
def FriedmanSTest(ResultsFunction, func_num, size, NRuns):
    rankarray = np.zeros((NRuns,size))    
    for i in range(NRuns):
        rankarray[i] = get_fract_ranks(np.transpose(ResultsFunction[:,i,func_num]))
    sumranks = np.zeros(size)
    avgranks = np.zeros(size)
    for i in range(size):
        sumranks[i] = np.sum(rankarray[:,i])
    avgranks = sumranks / NRuns
    raverage = (size+1)/2.0
    FriedmanS = 0
    for i in range(size):
        FriedmanS += (avgranks[i] - raverage)*(avgranks[i] - raverage)
    FriedmanS *= 12*NRuns/(size*(size+1.))
    #print(FriedmanS)
    return avgranks
def get_fract_ranks_and_groups(data):
    sort_index = np.argsort(-data)
    sort_list = -np.sort(-data)
    groups = []
    my_new_ranks = np.zeros(data.shape[0])
    counter = 0
    while(True):
        if(counter == data.shape[0]):
            break
        if(counter == data.shape[0]-1):
            my_new_ranks[counter] = counter
            break
        if(sort_list[counter] != sort_list[counter+1]):
            my_new_ranks[counter] = counter
            counter+=1            
        else:
            avgrank = 0
            start = counter
            while(sort_list[start] == sort_list[counter]):
                avgrank += counter
                counter+=1                
                if(counter == data.shape[0]):
                    break
            avgrank = avgrank / (counter - start)
            groups.append(counter - start)
            for i in range(start,counter):
                my_new_ranks[i] = avgrank
    index_rank = np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        new_rank_inv = data.shape[0] - my_new_ranks[i]
        index_rank[sort_index[i]] = new_rank_inv   
    return index_rank, groups
colors_s = [(0.0, 0.0, 0.0),(0.9, 0.0, 0.0),(0.0, 0.9, 0.0),(0.0, 0.0, 0.9),(0.9, 0.9, 0.0),(0.9, 0.0, 0.9),(0.0, 0.9, 0.9),(0.5, 0.5, 0.0),(0.0, 0.5, 0.5)]
markers = ['o','v','^','*','d','s','+','x','P']


PATH = ""
names = ["BlockEA","IEACOP","RDE","mLSHADE_LR","L_SRTDE","jSOa"]
AllRes = np.zeros((6,29,25,1001))
for func in range(29):    
    AllRes[0,func] = np.loadtxt(PATH+names[0]+"_"+str(func+1)+".txt").T
    AllRes[1,func,:,:1000] = np.loadtxt(PATH+names[1]+"_"+str(func+1)+".txt").T
    AllRes[2,func] = np.loadtxt(PATH+names[2]+"_"+str(func+1)+".txt").T
    AllRes[3,func,:,:1000] = np.loadtxt(PATH+names[3]+"_"+str(func+1)+".txt").T
    AllRes[4,func] = np.loadtxt(PATH+names[4]+"_"+str(func+1)+".txt").T
    AllRes[5,func,:,:1000] = np.loadtxt(PATH+names[5]+"_"+str(func+1)+".txt").T

#check when algorithms reach 1e-8 error, save these values to 1001-th cell in an array
for alg in range(6):
    for func in range(29):
        for run in range(25):
            savedFE = 0
            for it in range(1000):
                if(AllRes[alg,func,run,it] < 1e-8):
                    AllRes[alg,func,run,it] = 0 #reached 1e-9, set to zero
                    if(savedFE == 0):
                        #print(alg,func,run,AllRes[alg,func,run,1000],300*(it+1))
                        AllRes[alg,func,run,1000] = 300*(it+1)
                        savedFE = 1 #do not update further
            if(savedFE == 0): #never reached 1e-8
                #print("no",alg,func,run,AllRes[alg,func,run,1000],300000)
                AllRes[alg,func,run,1000] = 300000


#Mann-Whitney comparison
NFuncs = 29
NRuns = 25
MaxFEval = 300000
for alg1index in range(0,6):
    for alg2index in range(4,5):
        total = 0
        total2 = 0
        nplus = 0
        nminus = 0
        neq = 0
        for FNum in range(0,NFuncs):        
            avg1FE = np.reshape(AllRes[alg1index,FNum,:,-1],NRuns)#final NFE
            avg1ER = np.reshape(AllRes[alg1index,FNum,:,-2],NRuns)#final error
            avg2FE = np.reshape(AllRes[alg2index,FNum,:,-1],NRuns)#final NFE
            avg2ER = np.reshape(AllRes[alg2index,FNum,:,-2],NRuns)#final error
            avg1 = np.zeros(25)
            avg2 = np.zeros(25)          
            for i in range(25):
                if(avg1FE[i] == MaxFEval):
                    avg1[i] = avg1ER[i]
                else:
                    avg1[i] = avg1FE[i]*10e-9/MaxFEval
                if(avg2FE[i] == MaxFEval):
                    avg2[i] = avg2ER[i]
                else:
                    avg2[i] = avg2FE[i]*10e-9/MaxFEval
            #print(FNum+1,avg1,avg2)
            test = MannWhitneyU_myZ(avg1,avg2)
            total += test[0]
            if(np.isnan(test[1]) == False):
                total2 += test[1]        
            if(test[0] == 1):
                nplus += 1
            if(test[0] == -1):
                nminus += 1
            if(test[0] == 0):
                neq += 1
            print("func",FNum+1,test,np.mean(avg1ER),np.mean(avg2ER),np.mean(avg1FE),np.mean(avg2FE))
        print()
        print(names[alg1index],'vs',names[alg2index],"(",nplus,"+/",neq,"=/",nminus,"-) total:",total,", Z-score:",total2)
    print()
    print()
    
    
#Friedman
tocompare = [0,1,2,3,4,5]
NTests = len(tocompare)
NDim = 1
dimindex = 0
FRTotal = np.zeros((NDim,NTests))  
DIM = 30      
for FNum in range(0,NFuncs):   
    if(FNum == 1):
        continue
    ResFunction = np.zeros((NTests,NRuns,NFuncs))    
    FRFunc = np.zeros(NTests) 
    for i2,alg2 in enumerate(tocompare):     
        avg1FE = np.reshape(AllRes[alg2,FNum,:,-1],NRuns)    
        avg1ER = np.reshape(AllRes[alg2,FNum,:,-2],NRuns)      
        avg1 = np.copy(avg1ER)
        for i in range(NRuns):
            if(avg1FE[i] == MaxFEval):
                avg1[i] = avg1ER[i]
            else:
                avg1[i] = avg1FE[i]*10e-9/MaxFEval   
        ResFunction[i2,:,FNum] = avg1
    FRFunc = FriedmanSTest(ResFunction, FNum, NTests, NRuns)+1
    FRTotal[dimindex,:] += FRFunc
    #print(sum(FRFunc))
print("Friedman ranking:")
tmp2 = ""
for i2,alg2 in enumerate(tocompare):
    tmp2+= names[alg2]  
    #for dimindex in range(0,NDim):
    #    tmp2 +=  " & "+' %.2f'%(FRTotal[dimindex,i2])+''
    tmp2 += " & "+' %.2f'%(np.sum(FRTotal[:,i2]))+''
    tmp2 += "\\\\\n"
    #tmp2 += "\\hline\n"
print(tmp2)



#U-scores
tocompare = [0,1,2,3,4,5]
NTests = len(tocompare)
NDim = 1
FRTotal = np.zeros((NDim,NTests))  
NRuns = 25
for dimindex in range(0,1):    
    #DIM = DIMs[dimindex]     
    for FNum in range(0,NFuncs):
        if(FNum == 1):
            continue
        for i1,alg1 in enumerate(tocompare):
            for i2,alg2 in enumerate(tocompare):
                if(i1 < i2): 
                    #print(alg1,alg2)
                    avg1FE = np.reshape(AllRes[alg1,FNum,:,-1],NRuns)
                    avg1ER = np.reshape(AllRes[alg1,FNum,:,-2],NRuns)
                    avg2FE = np.reshape(AllRes[alg2,FNum,:,-1],NRuns)
                    avg2ER = np.reshape(AllRes[alg2,FNum,:,-2],NRuns)
                    avg1 = np.zeros(NRuns)
                    avg2 = np.zeros(NRuns)
                    #print(avg1ER)
                    #print(avg2ER)
                    #print(avg1FE)
                    #print(avg2FE)
                    for i in range(NRuns):
                        if(avg1FE[i] == MaxFEval):
                            avg1[i] = avg1ER[i]
                            #print(i,"ER")
                        else:
                            avg1[i] = avg1FE[i]*10e-9/MaxFEval
                            #avg1[i] = avg1ER[i]
                            #print(i,"FE+ER",avg1FE[i],MaxFEval[dimtest-1])
                        if(avg2FE[i] == MaxFEval):
                            avg2[i] = avg2ER[i]
                            #print(i,"ER")
                        else:
                            avg2[i] = avg2FE[i]*10e-9/MaxFEval       
                            #avg2[i] = avg2ER[i]
                            #print(i,"FE+ER",avg1FE[i],MaxFEval[dimtest-1])
                        #avg1 = np.copy(avg1ER)
                        
                    Sample1 = avg1
                    Sample2 = avg2
                    NewSample = -np.concatenate((Sample1,Sample2),axis=0)
                    NewRanks, Groups = get_fract_ranks_and_groups(NewSample)
                    SumRanks = 0
                    SumRanks2 = 0
                    for i in range(Sample1.shape[0]):
                        SumRanks += NewRanks[i]
                        SumRanks2 += NewRanks[Sample1.shape[0]+i]
                    U1 = SumRanks - Sample1.shape[0]*(Sample1.shape[0]+1.0)/2.0
                    U2 = SumRanks2 - Sample2.shape[0]*(Sample2.shape[0]+1.0)/2.0
                    #print(i1,i2)
                    #print(Sample1,Sample2)
                    #print(NewRanks)
                    #print(FNum,U1,U2,np.median(Sample1),np.median(Sample2))
                    FRTotal[dimindex,i1] += U1
                    FRTotal[dimindex,i2] += U2                    
                    #print(SumRanks,SumRanks2)
    #print(FRTotal[dimindex])
print("U-scores:")
tmp2 = ""
for i2,alg2 in enumerate(tocompare):
    tmp2+= names[alg2]  
    #for dimindex in range(0,2):
    #    tmp2 +=  " & "+' %.1f'%(FRTotal[dimindex,i2])+''
    tmp2 += " & "+' %.1f'%(np.sum(FRTotal[:,i2]))+''
    tmp2 += "\\\\\n"
    #tmp2 += "\\hline\n"
print(tmp2)



ticks_graph = []
DIM = 30
new_ticks = np.zeros(1000)
for k in range(1000):
    new_ticks[k] = DIM*10*(k+1)
    if(k%100*DIM == 0):
        ticks_graph.append(str(int(new_ticks[k])))
        
MaxFE = 300000

fig = plt.figure(figsize=(30, 30))
gs = gridspec.GridSpec(6,5)
gs.update(wspace=0.25, hspace=0.45) # set the spacing between axes. 
ax = []
for i in range(29):
    ax.append(plt.subplot(gs[i]))        
for FNum in range(0,29):       
    lines = []            
    ax[FNum].grid(True)  

    alg1 = names[0]
    cut1 = AllRes[0][FNum]
    avg1 = np.zeros(1000)
    for i in range(1000):
        avg1[i] = np.mean(cut1[:,i])
    line1, = ax[FNum].plot(new_ticks, avg1, color=colors_s[0],
                     lw=1, linestyle='-')
    lines.append(line1)

    alg2 = names[1]
    cut2 = AllRes[1][FNum]
    avg2 = np.zeros(1000)
    for i in range(1000):
        avg2[i] = np.mean(cut2[:,i])
    line2, = ax[FNum].plot(new_ticks, avg2, color=colors_s[1],
                     lw=1, linestyle='-')
    lines.append(line2) 

    alg3 = names[2]
    cut3 = AllRes[2][FNum]
    avg3 = np.zeros(1000)
    for i in range(1000):
        avg3[i] = np.mean(cut3[:,i])
    line3, = ax[FNum].plot(new_ticks, avg3, color=colors_s[2],
                     lw=1, linestyle='-')
    lines.append(line3) 

    alg4 = names[3]
    cut4 = AllRes[3][FNum]
    avg4 = np.zeros(1000)
    for i in range(1000):
        avg4[i] = np.mean(cut4[:,i])
    line4, = ax[FNum].plot(new_ticks, avg4, color=colors_s[3],
                     lw=1, linestyle='-')
    lines.append(line4) 

    alg5 = names[4]
    cut5 = AllRes[4][FNum]
    avg5 = np.zeros(1000)
    for i in range(1000):
        avg5[i] = np.mean(cut5[:,i])
    line5, = ax[FNum].plot(new_ticks, avg5, color=colors_s[4],
                     lw=1, linestyle='-')
    lines.append(line5) 

    alg6 = names[5]
    cut6 = AllRes[5][FNum]
    avg6 = np.zeros(1000)
    for i in range(1000):
        avg6[i] = np.mean(cut6[:,i])
    line6, = ax[FNum].plot(new_ticks, avg6, color=colors_s[5],
                     lw=1, linestyle='-')
    lines.append(line6) 

    ax[FNum].set_yscale('log')
    #ax[FNum].set_xscale('log')
    ax[FNum].autoscale(enable=True, axis='y', tight=None)      
    #ax[FNum].set_xticks(ticks=new_ticks)
    #ax[FNum].set_xticklabels(labels=ticks_graph)
    ax[FNum].tick_params(axis='x', rotation=0)
    ax[FNum].set_xlabel(r"$NFE$")
    ax[FNum].xaxis.set_label_coords(1.0, -0.2)
    ax[FNum].set_ylabel(r"$f(x)$", rotation=0)
    ax[FNum].yaxis.set_label_coords(-0.05, 1.03)
    ax[FNum].set_title(r"F"+str(FNum+1)+"")
    ax[FNum].set_xlim(0,MaxFE)

    if(FNum == 0):
        matplotlib.rcParams.update({'font.size': 8})
        ax[FNum].legend(lines, (alg1,alg2,alg3,alg4,alg5,alg6))
        matplotlib.rcParams.update({'font.size': 10})

#fig.savefig("comparison_1"+".eps",bbox_inches='tight')
fig.savefig("comparison_1"+".png",bbox_inches='tight')
#fig.savefig("comparison_1"+".svg",bbox_inches='tight')
fig.savefig("comparison_1"+".pdf",bbox_inches='tight')
plt.show()