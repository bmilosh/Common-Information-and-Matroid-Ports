##############################################################################################
############## CI Improved (using only sets that gave bounds in the prev. run) ###############
##############################################################################################

from itertools import permutations
from time import localtime, strftime, time

from gurobipy import *

import config
import re
from fibonew2 import (CI, AccStrCompatiblemnew, AK2exp, Init1m, Resol2m, bi, sb, ib,
                      Shannonm, SolSetExtractor2InfoCI, setgenerator)
from timing import endlog, log


def LambdaBoundsImproved(allsets):
    start = time()
    log("Start Program")

    errorfile = open('ErrorfilesCI.txt','a+')
    errorfile.write('++++++++++++++++++++++\n')
    errorfile.write('Files that give errors.\n')
    errorfile.write('++++++++++++++++++++++\n')
    allsets = setgenerator(config.destiny)
    counter = 1
    nonCI = 0
    oneCI = 0

    sfiles = open('checked_filesCIIm.txt','a+')
    nm = 'File listing checked files (and their improved CI bounds) from current run.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))

    for key in allsets:
        begin = time()
        log("Start run for {} CI (Improved)".format(key)) 
        sp = str(key).split()
        key1 = '{}_p{}'.format(sp[0],sp[2])
        
        fname = open('{}CIIm.txt'.format(key1),"w+")
        checker = 0
        try:
            fname11 = open('{}CI.txt'.format(key1))
            fname12 = open('{}CI.txt'.format(key1))
            fname13 = open('{}CI.txt'.format(key1))
        except:
            checker = 1
        if checker == 1:
            errorfile.write('File {} does not exist\n'.format(key1))
            continue

        llist = SolSetExtractor2InfoCI(fname11,fname13)  
        if llist == []:
            fname.write('{} is 1-CI.'.format(key))
            print('empty list')
            print('{}is 1-CI.'.format(key))
            fname.close()
            continue

        for line in fname12:
            line = line.strip()
            cifind = re.findall('^La.+is ([1.0-9]+)',line)
            if len(cifind) > 0:
                cival = float(cifind[0][:len(cifind[0])-1])
        
        hd = 'Solution file for {} CI (Improved)'.format(key)
        hd1 = '%' * len(hd)
        fname.write('{}\n'.format(hd1))
        fname.write('{}\n'.format(hd))
        fname.write('{}\n'.format(hd1))

        Astr = allsets[key]
        fdict = dict()
        flist = list()
        
        for i in range(len(llist)):
            for j in range(i+1,len(llist)):
                config.p = Model("gurotest")
                config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
                Init1m()
                Shannonm()
                AccStrCompatiblemnew(Astr)
                S = llist[i]
                T = llist[j]
                print(S,T)
                CI( bi(sb(S[0])), bi(sb(S[1])), 2**config.Part )
                CI( bi(sb(T[0])), bi(sb(T[1])), 2**(config.Part + 1) )
                Resol2m()
                if config.p.status == GRB.Status.OPTIMAL and config.p.objVal >= cival + 0.0001:
                    fdict['{0}-{1} and {2}-{3}'.format(S[0], S[1], T[0], T[1])] = config.p.objVal
                    for thing in range(1,config.Part):
                        fname.write('p{0} --> {1}\n'.format(thing,config.w[2**thing].X))
                    fname.write('Aux variable 1 --> {}\n'.format(config.w[2**config.Part].X))
                    fname.write('Aux variable 2 --> {}\n'.format(config.w[2**config.Part+1].X))
                    fname.write('Solution sets {0}-{1} and {2}-{3} have optimal value {4}.\n'.format(S[0], S[1], T[0], T[1],config.p.objVal))
        if fdict == {}:
            oneCI += 1
            Largest = "Nothing yet"
            smallest = "Nothing yet"
            print('Matroid port {} is non CI-improvable.'.format(key))
            sfiles.write('Matroid port {} is non CI-improvable.\n'.format(key))
        else:
            fname.write('\nThe results:\n')
            for keys in fdict:
                flist.append(fdict[keys])
                fname.write('{} : {}\n'.format(keys,fdict[keys]))
            Largest = max(flist)
            smallest = min(flist)
            nonCI += 1
            print('Matroid port {} is CI-improvable with a {} bound on lambda.'.format(key,Largest))
            sfiles.write('Matroid port {} is CI-improvable with a {} bound on lambda.\n'.format(key,Largest))
        fname.write('\nLargest is {}.'.format(Largest))
        fname.write('\nSmallest is {}.'.format(smallest))
        fname.write('\nSolution dictonary has {} items.'.format(len(fdict)))
        fname.close()  
        endlog(begin)
        if counter < len(allsets):
            difference = len(allsets)-counter
            if difference > 1:
                print('{0}done. {1} ports remaining. Moving to the next one...'.format(key,difference))
            else:
                print('{}done. One port left.'.format(key))
        else:
            print('*********************************************************')
            print('Last run made. Program concluded.')
            print('*********************************************************')
            sfiles.write('\n All {} ports checked.\n'.format(len(allsets)))
            if nonCI == 0:
                sfiles.write('All {} ports are non CI-improvable.\n'.format(oneCI))
            elif nonCI == 1 and nonCI != len(allsets):
                if oneCI == 1:
                    sfiles.write('There is one CI-improvable and {} non CI-improvable port here.\n'.format(oneCI))
                else:
                    sfiles.write('There is one CI-improvable and {} non CI-improvable ports here.\n'.format(oneCI))
            elif nonCI > 1 and nonCI < len(allsets):
                if oneCI == 1:
                    sfiles.write('There are {} CI-improvable ports, and {} non CI-improvable port here.\n'.format(nonCI,oneCI))
                else:
                    sfiles.write('There are {} CI-improvable ports, and {} non CI-improvable ports here.\n'.format(nonCI,oneCI))
            elif nonCI == len(allsets):
                sfiles.write('All {} ports are CI-improvable.\n'.format(nonCI))
        counter += 1    
    endlog(start)
    errorfile.close()

#LambdaBoundsImproved(config.destiny)