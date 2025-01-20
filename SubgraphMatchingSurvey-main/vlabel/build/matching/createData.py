
import numpy as np
algo=["_FT0","_FT1","_FT2","_FT3","_FT4","_FT5"]
algo=["DV","FT1_DV","FT2_DV","FT3_DV","FT4_DV","FT5_DV","LFTJ","RM","VEQ"]
timts=["1","10","50","100","1000"]
#timts=["10"]
timtsn=[1,10,50,100,1000]
#timtsn=[10]
algo=["DV","FT1_DV","FT2_DV","FT3_DV","FT4_DV","LFTJ","RM","VEQ"]
#algo=["_DVGM","FT1_DVGM","FT2_DVGM","FT3_DVGM","FT4_DVGM","_LFSK","_LFTJK"]
algo=["DLSBS"]

#timtsn=[1]
#FMPath="RE_PL_-1_patents_G16"
FMPath="performance_experiment/NSCORRLFTJ_LFTJ_-1_dblp_G_8"
FMPath2="performance_experiment/Skase_VEQ_-1_dblp_G_8"
FMPath="performance_experiment/ADIVERSITY_F0__LFTJ_-1_dblp_G_8"
FMPath="performance_experiment/DIV-BASIC_"

#FMPath="RE_PL_-1_patents_G"
pathsize=["8","16","24","32","64"]
#pathsize=["8"]
endf=".csv"
#TOPK=["10","50","100","500","1000","2500","5000","10000"]
TOPK=[]
avg_MC = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC1 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC2 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC3 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC4 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC5 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC6 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
avg_MC7 = [[0 for _ in range(len(timts))] for _ in range(len(algo))]
MC=[]
MC1=[]
MC2=[]
MC3=[]
MC4=[]
MC5=[]
MC6=[]
MC7=[]
AMC=[]
AMC1=[]
AMC2=[]
AMC3=[]
AMC4=[]
AMC5=[]
AMC6=[]
AMC7=[]
print(np.shape(pathsize))
TOPK=["10","50","100","500","1000","2500","5000","10000"]
TOPK_N=["10","50","100","500","1000","2500","5000","10000"]

 #with open(f'{folderall_}/{at}/{i}/fusbal_Tresults{t}.txt', 'r') as file:
for i in range(len(algo)):
#for i in range(len(pathsize)):
    for j in range(len(timts)):
        counter=0
        #with open(f'{FMPath}_{timts[j]}{algo[i]}{endf}', 'r') as file:
        with open(f'{FMPath}{algo[i]}_-1_youtube_G_64_{timts[j]}{endf}', 'r') as file:
            for line in file:
                values = line.split()
                if(counter<10000):
                    MC.append(float(values[4]))
                    counter=counter+1
                MC1.append(float(values[11]))
                MC2.append(float(values[12]))
                MC3.append(float(values[13]))
                MC4.append(float(values[14]))
                MC5.append(float(values[15]))
                MC6.append(float(values[16]))
                MC7.append(float(values[17]))
        avg_MC[i][j]=np.mean(MC)
        avg_MC1[i][j]=np.mean(MC1)
        avg_MC2[i][j]=np.mean(MC2)
        avg_MC3[i][j]=np.mean(MC3)
        avg_MC4[i][j]=np.mean(MC4)
        avg_MC5[i][j]=np.mean(MC5)
        avg_MC6[i][j]=np.mean(MC6)
        avg_MC7[i][j]=np.mean(MC7)
        #avg_UN[i][j]=np.mean(UN)
       # avg_TC[i][j]=np.mean(TC)
        #avg_TP[i][j]=np.mean(TP)
        if(len(MC)!=100 and len(MC)!=200):
            print("and so we pray",len(MC))
            print(algo[i],timts[j])
        MC.clear()
        MC1.clear()
        MC2.clear()
        MC3.clear()
        MC4.clear()
        MC5.clear()
        MC6.clear()
        MC7.clear()

print(FMPath)



#for i in range(len(algo)):
#for i in range(len(pathsize)):
#    for j in range(len(timts)):
        #with open(f'{FMPath}_{timts[j]}{algo[i]}{endf}', 'r') as file:
#        with open(f'{FMPath2}_{timts[j]}{endf}', 'r') as file:
#            for line in file:
#                values = line.split()
 ##               AMC.append(float(values[3]))
#print(MC)
#print(AMC)
#for i in range(400):
#    if(MC[i]!=AMC[i]):
#        print(i)

i=0
print("Time ALG ALG10 ALG50 ALG100 ALG250 ALG500 ALG750 ALG1000")
print("TOPK F2DVGM F4DVGM DVGMSQ LFSK LFTJK")
for j in range(len(timts)): 
    for t in range(1):   
    #for t in range(len(TOPK_N)):
        print(timts[j],end=' ')
        if(t==0):
            for i in range(len(algo)):
                print(avg_MC[i][j],end=' ')
        if(t==1):
            for i in range(len(algo)):
                print(avg_MC1[i][j],end=' ')       
        if(t==2):
            for i in range(len(algo)):
                print(avg_MC2[i][j],end=' ')
        if(t==3):
            for i in range(len(algo)):
                print(avg_MC3[i][j],end=' ')      
        if(t==4):
            for i in range(len(algo)):
                print(avg_MC4[i][j],end=' ')
        if(t==5):
            for i in range(len(algo)):
                print(avg_MC5[i][j],end=' ')       
        if(t==6):
            for i in range(len(algo)):
                print(avg_MC6[i][j],end=' ')
        if(t==7):
            for i in range(len(algo)):
                print(avg_MC7[i][j],end=' ')    
        print()