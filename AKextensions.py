from gurobipy import *
from itertools import combinations
from time import localtime, strftime, time

import config
from fibonew2 import (
    AK2exp, InitMatNew, MatroidCompatible, Resol2m, bi, bs, disjoint, rankfinder,
    ib, sb)
from timing import endlog, log


def CheckOneAK(mbases,gset,rnk):
    '''
    We check if the matroid is 1-AK.
    We also implement some optimizations that help reduce the number of flats to check.

    mbases (dictionary) containing matroids to be checked and their bases
    gset (list) is the ground set of the matroids
    rnk (int) is the rank of the matroids
    
    Note: observe that the matroids in a particular run need to be of the same size and rank.
    '''
    start = time()
    log("Start Program")

    checker = 1
    nonAK = 0
    oneAK = 0
    noAKmats = list()

    sfiles = open('runresultAK.txt','a+')
    nm = 'File listing checked files from polymatroid extension run using all sets (AK) with some optimizations.'
    nm1 = '%' * len(nm)
    sfiles.write('{}\n'.format(nm1))
    sfiles.write('{}\n'.format(nm))
    sfiles.write('{}\n'.format(nm1))

    for key in mbases:
        counter = 0
        begin = time()
        log('Start polymatroid extension check (AK) for {} using all sets with some optimizations.'.format(key))
        rankd,allranks = rankfinder(mbases[key],gset)  
        Ar = list()
        for i in range(0,len(allranks)-1,2): 
            if len(allranks[i]) < 2: continue
            #if Ar1[i+1] == rnk: continue  # not sure why this is here

            Ar.append(set([str(it3) for it3 in allranks[i]]))
        combs3 = combinations( [i for i in range(len(Ar))], 3)
        comb_hlder = list()

        ################################################
        ## We remove tuples (U,V,Z) where: 
        ## (i) UV has full rank
        ## (ii) U and V are subsets of Z
        ## (iii) UV is a subset of Z
        ## (iv) Z is the intersection of U and V
        ## (v) Z is a subset of UV
        ## (vi) UV and Z are a modular pair
        ## (vii) UV and Z have zero mutual information
        ################################################
        for combo in combs3:
            pre_comb_hlder = list()

            cmbs12 = Ar[combo[0]].union(Ar[combo[1]])
            excld = set([int(itm) for itm in cmbs12])
            ind = allranks.index(excld)
            rnk_excld = allranks[ind + 1]
            
            if rnk_excld == rnk: continue 
            
            if (Ar[combo[0]].issubset(Ar[combo[2]]) and Ar[combo[1]].issubset(Ar[combo[2]])) or cmbs12.issubset(Ar[combo[2]]): continue 
            if Ar[combo[2]]==Ar[combo[0]].intersection(Ar[combo[1]]) or cmbs12.issuperset(Ar[combo[2]]): continue

            #int_combo01 = [int(item) for item in cmbs12]
            set_combo01 = set( [int(item) for item in cmbs12] )
            index_combo01 = allranks.index(set_combo01)
            rnk_combo01 = allranks[index_combo01+1]

            #int_combo2 = [int(item) for item in Ar[combo[2]]]
            set_combo2 = set( [int(item) for item in Ar[combo[2]]] )
            index_combo2 = allranks.index(set_combo2)
            rnk_combo2 = allranks[index_combo2+1]

            combo_inters = cmbs12.intersection(Ar[combo[2]])
            #int_combointers = [int(item) for item in combo_inters]
            set_combointers = set( [int(item) for item in combo_inters] )
            index_combointers = allranks.index(set_combointers)
            rnk_combointers = allranks[index_combointers+1]

            combo_union = cmbs12.union(Ar[combo[2]])
            #int_combounion = [int(item) for item in combo_union]
            set_combounion = set( [int(item) for item in combo_union] )
            index_combounion = allranks.index(set_combounion)
            rnk_combounion = allranks[index_combounion+1]

            check_modularity = rnk_combo01 + rnk_combo2 - rnk_combounion - rnk_combointers
            mutual_info = rnk_combo01 + rnk_combo2 - rnk_combounion

            if check_modularity != 0 and mutual_info != 0: 
                pre_comb_hlder.append(Ar[combo[0]])
                pre_comb_hlder.append(Ar[combo[1]])
                pre_comb_hlder.append(Ar[combo[2]])
                comb_hlder.append(pre_comb_hlder)
        print('{} has {} 3-member working combinations.'.format(key,len(comb_hlder)))
        
        for i in range(len(comb_hlder)):
            combo1 = comb_hlder[i]
            J = combo1[0]
            K = combo1[1]
            L = combo1[2]
                                
            config.p = Model("gurotest")
            config.w = config.p.addVars(range(0,2**config.vrbls+1),name="w")
            InitMatNew()
            MatroidCompatible(mbases[key],gset)
            AK2exp(bi(sb(J)), bi(sb(K)), bi(sb(L)), 2**(config.Part))  
            Resol2m()
            if config.p.status == GRB.Status.OPTIMAL: continue
            print('{} is a non-AK matroid with violating sets {}, {} and {}.'.format(key,J,K,L))
            sfiles.write('{} is a non-AK matroid with violating sets {}, {} and {}.\n'.format(key,J,K,L))
            noAKmats.append(key)
            counter = 1
            break  ###### To find ALL combinations that break AK, suppress this line #####
        
        if counter == 0:
            oneAK += 1
            sfiles.write('{} is an AK matroid.\n'.format(key))
        else:
            nonAK += 1
        
        endlog(begin)
        if checker < len(mbases):
            difference = len(mbases)-checker
            if difference > 1:
                print('{0}done. {1} matroids remaining. Moving to the next one... \n'.format(key,difference))
            else:
                print('{}done. One matroid left.'.format(key))
        else:
            print('*********************************************************')
            print('Last run made. Program concluded.')
            print('*********************************************************')
            sfiles.write('\n All {} matroids checked.\n'.format(len(mbases)))
            if nonAK == 0:
                sfiles.write('All {} matroids are AK.\n'.format(oneAK))
            else:
                sfiles.write('non_AK_mats = {}\n'.format(noAKmats))
                if nonAK == 1 and nonAK != len(mbases):
                    if oneAK == 1:
                        sfiles.write('There is one non-AK and {} AK matroid here.\n'.format(oneAK))
                    else:
                        sfiles.write('There is one non-AK and {} AK matroids here.\n'.format(oneAK))
                elif nonAK > 1 and nonAK < len(mbases):
                    if oneAK == 1:
                        sfiles.write('There are {} non-AK matroids, and {} AK matroid here.\n'.format(nonAK,oneAK))
                    else:
                        sfiles.write('There are {} non-AK matroids, and {} AK matroids here.\n'.format(nonAK,oneAK))
                elif nonAK == len(mbases):
                    sfiles.write('All {} matroids are non-AK.\n'.format(nonAK))
        checker += 1
    endlog(start)
