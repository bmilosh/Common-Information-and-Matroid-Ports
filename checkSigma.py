#####################################################################
############## AK using all possible set combinations ###############
#####################################################################

from itertools import permutations
from time import localtime, strftime, time

from gurobipy import *

import config
from fibonew2 import (CI, AccStrCompatiblemnew, AK2exp, Init1m, Resol2m,
                      Shannonm, bs, disjoint, ib, setgenerator)
from timing import endlog, log


def SigmaBounds(onedict):
    start = time()
    log("Start Program")
    fru = setgenerator(onedict)

    nonAK = 0
    oneAK = 0
    noAKmats = list()
    counter = 1

    sfiles = open('runresult_sigma.txt','a+')
    nm = 'File listing bounds on sigma for current run.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))
    
    for key in fru:
        begin = time()
        log('Start run for {}AK (All sets).'.format(key))
        spl = str(key).split()
        
        fname = open('{0}_p{1}AK.txt'.format(spl[0],spl[2]),"w+")
        nm = 'Solution file for {} AK (All sets)\n'.format(key)
        nm1 = '%' * len(nm)
        fname.write('{}\n'.format(nm1))
        fname.write('{}\n'.format(nm))
        fname.write('{}\n'.format(nm1))

        fdict = dict()
        AStr = fru[key]  
        
        for i in range(1,2**config.Part+1):
            if len(bs(ib(i))) < 2: continue
            
            for j in range(i+1,2**config.Part+1):
                #if not disjoint(i,j): continue
                if len(bs(ib(j))) < 2: continue
                
                for k in range(j+1,2**config.Part+1):
                    #if not disjoint(i+j,k): continue
                    if len(bs(ib(k))) < 2: continue

                    config.p = Model("gurotest")
                    config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
                    Init1m()
                    Shannonm()
                    AccStrCompatiblemnew(AStr)
                    AK2exp(i,j,k,2**(config.Part))  
                    Resol2m()
                    if config.p.status == GRB.Status.OPTIMAL and config.p.objVal >= 1.01:
                        fdict['{0}, {1} and {2}'.format(bs(ib(i)),bs(ib(j)),bs(ib(k)))] = config.p.objVal
                        for thing in range(1,config.Part):
                            fname.write('p{0} --> {1}\n'.format(thing,config.w[2**thing].X))
                        fname.write('Aux variable --> {}\n'.format(config.w[2**config.Part].X))
                        fname.write('Solution sets {0}, {1} and {2} have optimal value {3}.\n'.format(bs(ib(i)),bs(ib(j)),bs(ib(k)),config.p.objVal))
        flist = list()
        if fdict == {}:
            Largest = "Nothing yet"
            smallest = "Nothing yet"
            print('Matroid port {} is 1-AK (bound on sigma is 1).'.format(key))
            sfiles.write('Matroid port {} is 1-AK (bound on sigma is 1).\n'.format(key))
            oneAK += 1
        else:
            fname.write('\nThe results:\n')
            for keys in fdict:
                flist.append(fdict[keys])
                fname.write('\'{0}\' : {1}\n'.format(keys,fdict[keys]))
            Largest = max(flist)
            smallest = min(flist)
            print('Matroid port {} has a {} bound on sigma.'.format(key,Largest))
            sfiles.write('Matroid port {} has a {} bound on sigma.\n'.format(key,Largest))
            noAKmats.append(key)
            nonAK += 1

        fname.write('\nLargest is {}.'.format(Largest))
        fname.write('\nSmallest is {}.'.format(smallest))
        fname.write('\nSolution dictonary has {} items.'.format(len(fdict)))
        fname.close()
        endlog(begin)
        
        if counter < len(onedict):
            difference = len(onedict)-counter
            if difference > 1:
                print('{0}done. {1} ports remaining. Moving to the next one... \n'.format(key,difference))
            else:
                print('{}done. One port left.'.format(key))
        else:
            print('*********************************************************')
            print('Last run made. Program concluded.')
            print('*********************************************************')
            sfiles.write('\n All {} matroid ports checked.\n'.format(len(onedict)))
            if nonAK == 0:
                sfiles.write('All {} matroid ports are AK.\n'.format(oneAK))
            else:
                sfiles.write('non_AK_mats = {}\n'.format(noAKmats))
                if nonAK == 1 and nonAK != len(onedict):
                    if oneAK == 1:
                        sfiles.write('There is one non-AK and {} AK matroid port here.\n'.format(oneAK))
                    else:
                        sfiles.write('There is one non-AK and {} AK matroid ports here.\n'.format(oneAK))
                elif nonAK > 1 and nonAK < len(onedict):
                    if oneAK == 1:
                        sfiles.write('There are {} non-AK matroid ports, and {} AK matroid port here.\n'.format(nonAK,oneAK))
                    else:
                        sfiles.write('There are {} non-AK matroid ports, and {} AK matroid ports here.\n'.format(nonAK,oneAK))
                elif nonAK == len(onedict):
                    sfiles.write('All {} matroid ports are non-AK.\n'.format(nonAK))
        counter += 1
    endlog(start)

#SigmaBounds(config.destiny)