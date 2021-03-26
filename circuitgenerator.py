# The tools we need to genrate the matroid ports of a given matroid are presented here.

from itertools import combinations


def mat_circ_gen(Q,B):
    '''
    We define a method to generate circuits of a matroid from its bases. 
    Algorithm appears in "Finding the circuits of a matroid" by E. Minieka, Journal of Research of the National Bureau of Standards (1976).

    # The idea is the following: For every basis B and every element x not in B, set C = {x}. 
    # If B\cup x\setminus x_i, for i\in\{1,...,r\}, is a basis, update C\cup x_i. 
    # C is a circuit of the matroid M.
   
    Q (list) is the ground set of the matroid.
    B (list) is the matroid's bases family.

    output: allcir (list) circuits of the matroid.
    '''
    allcir = list()
    frbases = list()
    for bas in B:
        #sbas = sorted(bas)
        #jbas = ''.join(sbas)
        frbases.append(''.join(sorted(bas)))
    for basis in frbases:
        for x in Q:
            cir = list()
            if basis.__contains__(x): continue
            cir.append(x)
            for i in range(len(basis)):
                #newbasis = basis.replace(basis[i],x)
                #soneba = sorted(newbasis)
                jonebas = ''.join(sorted( basis.replace(basis[i],x) ))
                if jonebas not in frbases: continue
                cir.append(basis[i])
            if len(cir) <= 1: continue
            if sorted(cir) in allcir: continue
            allcir.append(sorted(cir))
            #allcir.append(''.join( sorted(cir) ))
    allcircs = list()
    for thing in allcir:
        circuit = ''.join(thing)
        allcircs.append(circuit)
    return allcircs  


# A method to create min Gamma from the port of a matroid at a point.
def acc_str_from_matr_port(Q,B,x):
    '''
    A circuit forms a minimal authorised set of Gamma if that circuit contains the port x.
    Q (list) is the ground set of the matroid.
    B (list) is the matroid's bases family.
    x (int) is the given port.

    output: BB (list) minimal authorised sets of Gamma.
    '''
    dre = mat_circ_gen(Q,B)
    BB = list()
    for circuit in dre:
        if x not in circuit: continue
        bb = circuit.replace(x,'')
        BB.append(bb)
    return BB


def ASgenfrommatroidport(groundset,basfam): 
    '''
    This generates all the matroid ports corresponding to a matroid M.

    groundset (list) is the ground set of the matroids.
    basfam (dictionary) is a dictionary where each key-value pair is a matroid-bases family pair.
    Note that all the matroids should have the same ground set.
    '''
    groundset = [str(x) for x in groundset]
    allaccstrs = dict()
    for key in basfam:
        basesfamily = basfam[key]
        naccessstructure = list()
        for x in groundset:
            port = x
            accessstructure = acc_str_from_matr_port(groundset,basesfamily,port)
            naccessstructure_1 = list()
            if x == 0: continue
            for item in accessstructure:
                newitem = item.replace('0',x)
                naccessstructure_1.append(newitem)
            naccessstructure.append(naccessstructure_1)
            naccessstructure_1.sort()
            allaccstrs['{0} port {1}'.format(key,x)] = naccessstructure_1
    chker = 0
    print('acc_strs = {')
    for item in allaccstrs:
        
        if chker < len(allaccstrs) - 1:
            print('\'{}\' : {} ,'.format(item,allaccstrs[item]))  
            chker += 1
        else:
            print('\'{}\' : {} '.format(item,allaccstrs[item]))  
    print('}')
    return None 
