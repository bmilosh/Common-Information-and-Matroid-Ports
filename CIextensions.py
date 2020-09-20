from gurobipy import *
from itertools import combinations
from time import localtime, strftime, time

import config
from fibonew2 import (
    CI, AK2exp, InitMatNew, MatroidCompatible, Resol2m, bi, bs, disjoint, ClosedSetExtractor, rankfinder,
    ib, sb)
from timing import endlog, log


def CheckOneCI(mbases,gset):  
    '''
    We check if the matroid is 1-CI.
    This check is done using flats of the matroid.
    We also implement some optimizations that help reduce the number of flats to check.

    mbases (dictionary) containing matroids to be checked and their bases
    gset (list) is the ground set of the matroids

    Note: observe that the matroids in a particular run need to be of the same size.
    '''
    start = time()
    log("Start Program")

    checker = 1
    nonCI = 0
    oneCI = 0
    noCImats = list()
    
    sfiles = open('runresultCI.txt','a+')
    nm = 'File listing checked files from polymatroid extension run using all closed sets.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))

    for key in mbases:
        counter = 0 
        begin = time()
        log('Start polymatroid extension check for {} using all closed sets (with optimization).'.format(key))
        smrnks,allranks = rankfinder(mbases[key],gset) 
        Ar1 = ClosedSetExtractor(mbases[key],gset)   
        Ar = list()
        for item in Ar1:
            if len(item) < 2: continue  
            Ar.append( set([str(it3) for it3 in item]) )
        
        lAr = [i for i in range(len(Ar))]
        combs2 = combinations(lAr,2)
        comb_hlder = list()

        ##########################################################################
        ## We remove modular flats and also flats with zero mutual information. ##
        ##########################################################################
        for combo in combs2:
            pre_comb_hlder = list()

            set_combo0 = set( [int(item) for item in Ar[combo[0]]] )
            index_combo0 = allranks.index(set_combo0)
            rnk_combo0 = allranks[index_combo0+1]

            set_combo1 = set( [int(item) for item in Ar[combo[1]]] )
            index_combo1 = allranks.index(set_combo1)
            rnk_combo1 = allranks[index_combo1+1]

            combo_inters = Ar[combo[0]].intersection(Ar[combo[1]])
            set_combointers = set( [int(item) for item in combo_inters] )
            index_combointers = allranks.index(set_combointers)
            rnk_combointers = allranks[index_combointers+1]

            combo_union = Ar[combo[0]].union(Ar[combo[1]])
            set_combounion = set( [int(item) for item in combo_union] )
            index_combounion = allranks.index(set_combounion)
            rnk_combounion = allranks[index_combounion+1]

            check_modularity = rnk_combo0 + rnk_combo1 - rnk_combounion - rnk_combointers

            mutual_information = rnk_combo0 + rnk_combo1 - rnk_combounion

            if Ar[combo[0]].isdisjoint(Ar[combo[1]]) and check_modularity != 0 and mutual_information != 0:
                pre_comb_hlder.append(Ar[combo[0]])
                pre_comb_hlder.append(Ar[combo[1]])
                comb_hlder.append(pre_comb_hlder)
        print('{} has {} closed sets and {} 2-member working combinations.'.format(key,len(Ar),len(comb_hlder)))

        for i in range(len(comb_hlder)):
            combo1 = comb_hlder[i]
            J = combo1[0]
            K = combo1[1]
            
            config.p = Model("gurotest")
            config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
            InitMatNew()
            MatroidCompatible(mbases[key],gset)
            CI(bi(sb(J)), bi(sb(K)), 2**(config.Part))  
            Resol2m()
            if config.p.status == GRB.Status.OPTIMAL: continue
            print('{} is a non-CI matroid with violating sets {} and {}.'.format(key,J,K))
            sfiles.write('{} is a non-CI matroid with violating sets {} and {}.\n'.format(key,J,K))
            noCImats.append(key)
            counter = 1
            break  ###### To find ALL combinations that break CI, suppress this line #####

        if counter == 0:
            oneCI += 1
            sfiles.write('{} is a CI matroid.\n'.format(key))
        else:
            nonCI += 1
        
        endlog(begin)
        if checker < len(mbases):
            difference = len(mbases)-checker
            if difference > 1:
                print('{0} done. {1} matroids remaining. Moving to the next one... \n'.format(key,difference))
            else:
                print('{} done. One matroid left.'.format(key))
        else:
            print('*********************************************************')
            print('Last run made. Program concluded.')
            print('*********************************************************')
            sfiles.write('\n All {} matroids checked.\n'.format(len(mbases)))
            if nonCI == 0:
                sfiles.write('All {} matroids are CI.\n'.format(oneCI))
            else:
                sfiles.write('non_CI_mats = {}\n'.format(noCImats))
                if nonCI == 1 and nonCI != len(mbases):
                    if oneCI == 1:
                        sfiles.write('There is one non-CI and {} CI matroid here.\n'.format(oneCI))
                    else:
                        sfiles.write('There is one non-CI and {} CI matroids here.\n'.format(oneCI))
                elif nonCI > 1 and nonCI < len(mbases):
                    if oneCI == 1:
                        sfiles.write('There are {} non-CI matroids, and {} CI matroid here.\n'.format(nonCI,oneCI))
                    else:
                        sfiles.write('There are {} non-CI matroids, and {} CI matroids here.\n'.format(nonCI,oneCI))
                elif nonCI == len(mbases):
                    sfiles.write('All {} matroids are non-CI.\n'.format(nonCI))
        checker += 1
    endlog(start)
