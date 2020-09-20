##############################################################################################
############## AK Improved (using only sets that gave bounds in the prev. run) ###############
##############################################################################################

import re
from itertools import combinations, permutations
from time import localtime, strftime, time

from gurobipy import *

import config
from fibonew2 import (CI, AccStrCompatiblemnew, AK2exp, Init1m, Resol2m,
                      Shannonm, SolSetExtractor2InfoAK, bi, bs, ib, sb,
                      setgenerator)
from timing import endlog, log


def SigmaBoundsImproved(onedict):
    start = time()
    log("Start Program")

    errorfile = open('ErrorfilesAK.txt','a+')
    errorfile.write('++++++++++++++++++++++\n')
    errorfile.write('Files that give errors.\n')
    errorfile.write('++++++++++++++++++++++\n')

    fru = setgenerator(onedict)
    counter = 1
    nonAK = 0
    oneAK = 0

    sfiles = open('checked_filesAKIm.txt','a+')
    nm = 'File listing checked files (and their improved AK bounds) from current run.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))
    for key in fru:
        begin = time()
        log('Start run for {}AK (Improved).'.format(key))
        spl = str(key).split()
        key1 = '{}_p{}'.format(spl[0],spl[2])
        fname = open('{0}AKIm.txt'.format(key1),'w+')

        checker = 0
        try:
            sol_file = open('{}AK.txt'.format(key1))
            sol_file1 = open('{}AK.txt'.format(key1))
            sol_file2 = open('{}AK.txt'.format(key1))
        except:
            checker = 1
        if checker == 1:
            errorfile.write('File {} does not exist\n'.format(key1))
            continue

        llist = SolSetExtractor2InfoAK(sol_file,sol_file2)
        if llist == []:
            fname.write('{}is 1-AK.'.format(key))
            print('empty list')
            print('{}is 1-AK.'.format(key))
            fname.close()
            continue
        
        for line in sol_file1:
            line = line.strip()
            cifind = re.findall('^La.+is ([1.0-9]+)',line)
            if len(cifind) > 0:
                cival = float(cifind[0][:len(cifind[0])-1])

        hd = 'Solution file for {} AK (Improved)'.format(key)
        hd1 = '%' * len(hd)
        fname.write('{}\n'.format(hd1))
        fname.write('{}\n'.format(hd))
        fname.write('{}\n'.format(hd1))

        fdict = dict()
        AStr = fru[key]  
        for i in range(0,len(llist),3): 
            for j in range(i+3,len(llist),3): 
                if j==i: continue
                config.p = Model("gurotest")
                config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
                Init1m()
                Shannonm()
                AccStrCompatiblemnew(AStr)
                S11 = llist[i]
                S21 = llist[i+1]
                S31 = llist[i+2]
                T11 = llist[j]
                T21 = llist[j+1]
                T31 = llist[j+2]
                S12 = bi(sb(S11))
                S22 = bi(sb(S21))
                S32 = bi(sb(S31))
                T12 = bi(sb(T11))
                T22 = bi(sb(T21))
                T32 = bi(sb(T31))
                AK2exp(S12,S22,S32,2**(config.Part))
                AK2exp(T12,T22,T32,2**(config.Part + 1))  
                Resol2m()
                if config.p.status == GRB.Status.OPTIMAL and config.p.objVal >= cival + 0.0001:
                    fdict['{0}-{1}-{2} and {3}-{4}-{5}'.format(S11,S21,S31,T11,T21,T31)] = config.p.objVal
                    for thing in range(1,config.Part):
                        fname.write('p{0} --> {1}\n'.format(thing,config.w[2**thing].X))
                    fname.write('Aux variable 1 --> {}\n'.format(config.w[2**config.Part].X))
                    fname.write('Aux variable 2 --> {}\n'.format(config.w[2**config.Part+1].X))
                    fname.write('Solution sets {0}-{1}-{2} and {3}-{4}-{5} have optimal value {6}.\n'.format(S11,S21,S31,T11,T21,T31,config.p.objVal))
        flist = list()
        if fdict == {}:
            oneAK += 1
            Largest = "Nothing yet"
            smallest = "Nothing yet"
            print('Matroid port {} is non AK-improvable.'.format(key))
            sfiles.write('Matroid port {} is non AK-improvable.\n'.format(key))
        else:
            fname.write('\nThe results:\n')
            for keys in fdict:
                flist.append(fdict[keys])
                fname.write('\'{0}\' : {1}\n'.format(keys,fdict[keys]))
            Largest = max(flist)
            smallest = min(flist)
            nonAK += 1
            print('Matroid port {} is AK-improvable with a {} bound on sigma.'.format(key,Largest))
            sfiles.write('Matroid port {} is AK-improvable with a {} bound on sigma.\n'.format(key,Largest))
        fname.write('\nLargest is {}.'.format(Largest))
        fname.write('\nSmallest is {}.'.format(smallest))
        fname.write('\nSolution dictonary has {} items.'.format(len(fdict)))
        fname.close()
        endlog(begin)
        if counter < len(onedict):
            difference = len(onedict)-counter
            if difference > 1:
                print('{0}done. {1} ports remaining. Moving to the next one...'.format(key,difference))
            else:
                print('{}done. One port left.'.format(key))
        else:
            print('*********************************************************')
            print('Last run made. Program concluded.')
            print('*********************************************************')
            sfiles.write('\n All {} ports checked.\n'.format(len(onedict)))
            if nonAK == 0:
                sfiles.write('All {} ports are non AK-improvable.\n'.format(oneAK))
            elif nonAK == 1 and nonAK != len(onedict):
                if oneAK == 1:
                    sfiles.write('There is one AK-improvable and {} non AK-improvable port here.\n'.format(oneAK))
                else:
                    sfiles.write('There is one AK-improvable and {} non AK-improvable ports here.\n'.format(oneAK))
            elif nonAK > 1 and nonAK < len(onedict):
                if oneAK == 1:
                    sfiles.write('There are {} AK-improvable ports, and {} non AK-improvable port here.\n'.format(nonAK,oneAK))
                else:
                    sfiles.write('There are {} AK-improvable ports, and {} non AK-improvable ports here.\n'.format(nonAK,oneAK))
            elif nonAK == len(onedict):
                sfiles.write('All {} ports are AK-improvable.\n'.format(nonAK))
        counter += 1
    endlog(start)
    errorfile.close()

#SigmaBoundsImproved(config.destiny)