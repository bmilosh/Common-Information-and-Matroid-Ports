#####################################################################
############## CI using all possible set combinations ###############
#####################################################################

from itertools import permutations
from time import localtime, strftime, time

from gurobipy import *

import config
from fibonew2 import (CI, AccStrCompatiblemnew, AK2exp, Init1m, Resol2m,
                      Shannonm, bs, disjoint, ib, setgenerator)
from timing import endlog, log


def LambdaBounds(onedict):
    start = time()
    log("Start Program")

    fru = setgenerator(onedict)
    nonCI = 0
    oneCI = 0
    noCImats = list()
    counter = 1

    sfiles = open('runresult_lambda.txt','a+')
    nm = 'File listing bounds on lambda for current run.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))

    for key in fru:
        key1 = str(key).strip()
        spl = key1.split()
        begin = time()
        log('Start run for {} CI (All sets).'.format(key1))

        fname = open('{0}_p{1}CI.txt'.format(spl[0],spl[2]),"w+")
        nm = 'Solution file for {} CI (All sets)\n'.format(key)
        nm1 = '%' * len(nm)
        fname.write('{}\n'.format(nm1))
        fname.write('{}\n'.format(nm))
        fname.write('{}\n'.format(nm1))

        fdict = dict()
        AStr = fru[key]  
        
        for i in range(1,2**config.Part+1):
            for j in range(i+1,2**config.Part+1):
                if not disjoint(i,j): continue
                if len(bs(ib(i))) < 2 or len(bs(ib(j))) < 2: continue
                config.p = Model("gurotest")
                config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
                Init1m()
                Shannonm()
                AccStrCompatiblemnew(AStr)
                CI(i,j,2**config.Part)  
                Resol2m()
                if config.p.status == GRB.Status.OPTIMAL and config.p.objVal >= 1.01:
                    fdict['{0} and {1}'.format(bs(ib(i)),bs(ib(j)))] = config.p.objVal
                    for thing in range(1,config.Part):
                        fname.write('p{0} --> {1}\n'.format(thing,config.w[2**thing].X))
                    fname.write('Aux variable --> {}\n'.format(config.w[2**config.Part].X))
                    fname.write('Solution sets {0} and {1} have optimal value {2}.\n'.format(bs(ib(i)),bs(ib(j)),config.p.objVal))
        
        flist = list()
        if fdict == {}:
            Largest = "Nothing yet"
            smallest = "Nothing yet"
            print('Matroid port {} is 1-CI (bound on lambda is 1).'.format(key))
            sfiles.write('Matroid port {} is 1-CI (bound on lambda is 1).\n'.format(key))
            oneCI += 1
        else:
            fname.write('\nThe results:\n')
            for keys in fdict:
                flist.append(fdict[keys])
                fname.write('\'{0}\' : {1}\n'.format(keys,fdict[keys]))
            Largest = max(flist)
            smallest = min(flist)
            print('Matroid port {} has a {} bound on lambda.'.format(key,Largest))
            sfiles.write('Matroid port {} has a {} bound on lambda.\n'.format(key,Largest))
            noCImats.append(key)
            nonCI += 1
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
            if nonCI == 0:
                sfiles.write('All {} matroid ports are CI.\n'.format(oneCI))
            else:
                sfiles.write('non_CI_mats = {}\n'.format(noCImats))
                if nonCI == 1 and nonCI != len(onedict):
                    if oneCI == 1:
                        sfiles.write('There is one non-CI and {} CI matroid port here.\n'.format(oneCI))
                    else:
                        sfiles.write('There is one non-CI and {} CI matroid ports here.\n'.format(oneCI))
                elif nonCI > 1 and nonCI < len(onedict):
                    if oneCI == 1:
                        sfiles.write('There are {} non-CI matroid ports, and {} CI matroid port here.\n'.format(nonCI,oneCI))
                    else:
                        sfiles.write('There are {} non-CI matroid ports, and {} CI matroid ports here.\n'.format(nonCI,oneCI))
                elif nonCI == len(onedict):
                    sfiles.write('All {} matroid ports are non-CI.\n'.format(nonCI))
        counter += 1
    endlog(start)

#LambdaBounds(config.destiny)
