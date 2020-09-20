from gurobipy import *
import config
import re
from itertools import combinations, permutations

# bs: Converts a binary vector into a set.
def bs(v):
    '''
    v (list)

    output: s (set) set representation of v
    '''
    s=set([])
    for j in range(0,config.vrbls):
        if v[j]==1:
            s.add(j)
    return s

# bi: Converts a binary vector into an integer.
def bi(v):
    '''
    v (list)

    output: i (list) integer representation of v
    '''
    i=0
    for j in range(0,config.vrbls):
        i+=v[j]*2**j
    return i

# ib: Converts an integer into a binary vector.
def ib(n):
    '''
    n (int)

    output: v (list) binary representation of n
    '''
    v=[]
    i=1
    for i in range(1,config.vrbls+1):
        v=v+[n%2]
        n=(n-n%2)/2
        i+=1
    return v

# sb: Converts a set into a binary vector.
def sb(s):
    '''
    s (set or list)

    output: u (list) binary representation of s
    '''
    s = [int(it) for it in s]
    u=[0]*config.vrbls
    for i in range(0,config.vrbls): 
    	if i in s:
    		u[i]=1
    return u

# union: Similar to the union method of set, outputs a set union given two integers.
def union(i,j):
    '''
    i (int)
    j (int)

    output: (set) union of i and j
    '''
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls):
        if v[k]==1:
            u[k]=1
    return(bi(u))

# contained: Similar to the subset method of set, given two integers, checks if one is a subset of the other and returns a boolean.
def contained(i,j):
    '''
    i (int)
    j (int)

    output: boolean --> True if i contained in j; False otherwise
    '''
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls):
        if u[k]==1 and v[k]==0:
            return False
    return True

# disjoint: Similar to the disjoint method of set, given two integers, checks if both are disjoint and returns a boolean.
def disjoint(i,j):
    '''
    i (int)
    j (int)

    output: boolean --> True if disjoint; False otherwise
    '''
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls): # [0..vars-1]:
        if u[k]==v[k] and u[k]==1:
            return False
    return True

# Init1m: To initialise when searching for bounds.
def Init1m():
    config.p.setObjective(config.w[2**config.vrbls],GRB.MINIMIZE)
    for i in range(1,config.Part): 
    	config.p.addConstr(config.w[2**config.vrbls]>=config.w[2**i])
    config.p.addConstr(config.w[0]==0)
    config.p.addConstr(config.w[1]==1)

# Resol2m: Calls the optimize method to solve the LP.
def Resol2m():
    config.p.setParam('OutputFlag',0) # set to 0 to suppress the printing of the optimization; set to 1 otherwise
    res=config.p.optimize()
    return res

# rankfinder: Computes the ranks of all subsets of the matroid.
def rankfinder(matbases,groundset):
    '''
    We compute the ranks of all the subsets of a matroid's ground set
    matbases (list) the bases family of the matroid
    groundset (list) ground set of the matroid (each element is an integer)

    outputs: rankdict (dictionary) set-rank pair
             ranklist (list) a set followed by its rank
    '''
    matrank = len(matbases[0])
    subsetstuple = sum([list(combinations(groundset, i)) for i in range(len(groundset)+1)], [])
    allsubsets = [set(i) for i in subsetstuple]
    setlist = list() 
    rankdict = dict()
    ranklist = list()
    ##################################################
    ### We compute the rank of all independent subsets
    ##################################################
    for item in matbases:
        setlist1 = sum([list(combinations(item, i)) for i in range(len(item)+1)], [])
        for thing in setlist1:
            tuptoset = set([int(thi) for thi in thing])
            if tuptoset in setlist: continue
            rank = len(tuptoset)
            setlist.append(tuptoset)
            rankdict['{}'.format(tuptoset)] = rank
            ranklist.append(tuptoset)
            ranklist.append(rank)
    ##################################################
    ### We compute the rank of all dependent subsets
    ##################################################
    
    for asubset in allsubsets:
        if asubset in setlist: continue
        gh = list()
        chk = 0
        for items in setlist:
            gh2 = asubset.intersection(items)
            gh1 = len(gh2)
            if gh1 == matrank:
                chk = 1
                rankdict['{}'.format(asubset)] = gh1
                ranklist.append(asubset)
                ranklist.append(gh1)
                break
            gh.append(gh1)
        if chk == 1: continue
        de = max(gh)
        rankdict['{}'.format(asubset)] = de
        ranklist.append(asubset)
        ranklist.append(de)
    return rankdict, ranklist  

# MatroidCompatible: Ensures that we are dealing with a matroid (when searching for extensions).
def MatroidCompatible(Mbases,groundset):
    '''
    mbases (list) the bases family of the matroid
    groundset (list) ground set of the matroid (each element is an integer)
    '''
    subsetstuple = sum([list(combinations(groundset, i)) for i in range(len(groundset)+1)], [])
    subsetslist = [set(i) for i in subsetstuple]  # list()
    ranks,rankl = rankfinder(Mbases,groundset)
    for A in subsetslist:
        j=sum([2**i for i in A])
        config.p.addConstr(config.w[j] == ranks['{}'.format(set(A))])  

# ClosedSetExtractor: Generates all the closed sets of given matroid.
def ClosedSetExtractor(mbases,groundset):
    '''
    We extract the closed sets of the given matroid.
    A set is closed if the rank of each of its supersets is greater than its own rank.
    mbases (list) the bases family of the matroid
    groundset (list) ground set of the matroid (each element is an integer)

    output: finlist (list) all closed sets of the matroid
    '''
    ranks,rankl = rankfinder(mbases,groundset)
    lid = [key for key in ranks]
    finlst = list()
    
    for i in range(len(lid)):
        
        itl = set([int(item) for item in lid[i] if item not in '}{( )\',setfrozen']) 
        chckr = 0
        for j in range(i+1,len(lid)):
            jtl = set([int(item) for item in lid[j] if item not in '}{ ()\',setfrozen']) 
            if itl.isdisjoint(jtl):  
                continue
            if itl.intersection(jtl) != itl:   
                continue
            if ranks[lid[j]] == ranks[lid[i]]: 
                chckr = 0
                break
            chckr += 1
        if chckr > 0:
            finlst.append(itl)
    return finlst

# InitMatNew: To initialise when searching for extensions.
def InitMatNew():
    config.p.setObjective(0,GRB.MINIMIZE)	
    for i in range(0,config.vrbls):  
        config.p.addConstr(config.w[2**i]>=0)  # entropy of each singleton is non-negative
    for i in range(0,config.vrbls):  
        config.p.addConstr(config.w[2**config.vrbls-1]>=config.w[2**config.vrbls-1-2**i])  # entropy of groundset > groundset less an element
    #############################
    ## Submodularity condition ##
    #############################
    for i in range(0,2**config.vrbls):   
        v=ib(i)
        for j in range(0,config.vrbls):  
            if v[j]==0:
                for k in range(j+1,config.vrbls):  
                    if v[k]==0:
                        config.p.addConstr(config.w[i+2**j]+config.w[i+2**k]>=config.w[i]+config.w[i+2**j+2**k])  # submodularity condition      


# Shannon0: H(S_i)>=0
def Shannon0():
    for i in range(0,config.vrbls):
        config.p.addConstr(config.w[2**i]>=0)

# Shannon1: H(S_Q)\geq H(S_{Q-i}) for all i
def Shannon1():
    for i in range(0,config.vrbls):
        config.p.addConstr(config.w[2**(config.vrbls)-1]>=config.w[2**(config.vrbls)-1-2**i])

# Shannon2: Implements the submodularty condition
def Shannon2():
    for i in range(0,2**config.vrbls):
        v=ib(i)
        for j in range(0,config.vrbls):
            if v[j]==0:
                for k in range(j+1,config.vrbls):
                    if v[k]==0:
                        config.p.addConstr(config.w[i+2**j]+config.w[i+2**k]>=config.w[i]+config.w[i+2**j+2**k])

def Shannonm():
	Shannon0()
	Shannon1()
	Shannon2()

# CE: The conditional entropy function.
def CE(v,i1,i2):
	i12=union(i1,i2)
	s=v[i12]-v[i2]
	return s

# MIC: The conditional mutual information function.
def MIC(v,i1,i2,i3):
	i13=union(i1,i3)
	i23=union(i2,i3)
	i123=union(i1,i23)
	s=v[i13]+v[i23]-v[i123]-v[i3]
	return s

# AK2exp: Implements the AK lemma.
def AK2exp(i,j,k,l): 
	config.p.addConstr(CE(config.w,l,union(i,j))==0)
	config.p.addConstr(CE(config.w,i,l)==CE(config.w,i,k))
	config.p.addConstr(CE(config.w,j,l)==CE(config.w,j,k))
	config.p.addConstr(CE(config.w,union(i,j),l)==CE(config.w,union(i,j),k))

# MI: The mutual information function.
def MI(v,i1,i2):
	i12=union(i1,i2)
	s=v[i1]+v[i2]-v[i12]
	return s

# CI: Implements the CI property.
def CI(i,j,k):
    config.p.addConstr(CE(config.w,k,i)==0)
    config.p.addConstr(CE(config.w,k,j)==0)
    config.p.addConstr(MI(config.w,i,j)==config.w[k])

def setgenerator(structure):
    '''
    structure is the dictionary containing the matroid ports
    '''
    onedictionary = dict()
    for key in structure:
        crowd = list()
        for item in structure[key]:
            h = 0
            for thing in item:
                h += 2**int(thing)
            crowd.append(int(h))
        onedictionary[key] = set(crowd)
    return onedictionary

# AccStrCompatiblemnew: Ensures that we are dealing with a matroid (when searching for bounds).
def AccStrCompatiblemnew(Ad):
    '''
    Generally, Ad is the output of the setgenerator function
    '''
    for j in range(0,2**config.Part-1,2): 
        t=1
        for k in Ad:
            if contained(k,j):
                t=0
        config.p.addConstr(config.w[j+1]-config.w[j]==t)

# SolSetExtractor2InfoCI: Extracts the set combinations that break CI for a particular matroid.
def SolSetExtractor2InfoCI(file1,file2):
    '''
    We use this to extract the combinations that break CI for a particular matroid for use in search for even better bounds.

    file1, file2 (file) two .txt files that contain the sets that gave bounds in a previous run. Both files are the same.
    '''
    emlist = list()
    for line in file1:
        line = line.strip()
        # Solution sets {0, 2} and {1, 3} have optimal value 1.3333333333333335.
        se1 = re.findall('^S.+ sets ([{][0-9, ]+[}])',line)
        de1 = re.findall('^S.+ and ([{][0-9, ]+[}])',line)
        if len(se1) > 0 and len(de1) > 0:
            selist = list()
            delist = list()
            for i in range(1,len(se1[0]),3):
                fe = int(se1[0][i])
                selist.append(fe)
            for i in range(1,len(de1[0]),3):
                fe1 = int(de1[0][i])
                delist.append(fe1)
            emlist.append([selist,delist])

    if emlist == []:
        for line in file2:
            line = line.strip()
            se = re.findall('^Solution sets set\(([\\[][0-9, ]+[]])',line)  
            de = re.findall('^Solution.+\(([\\[][0-9, ]+[]])',line)
            #Solution sets set([0, 2]) and set([1, 5]) have optimal value 1.33333333333.
            if len(se) > 0 and len(de) > 0:
                selist = list()
                delist = list()
                for i in range(1,len(se[0]),3):
                    fe = int(se[0][i])
                    selist.append(fe)
                for i in range(1,len(de[0]),3):
                    fe1 = int(de[0][i])
                    delist.append(fe1)
                emlist.append([selist,delist])
    return emlist

# SolSetExtractor2InfoAK: Extracts the set combinations that break AK for a particular matroid.
def SolSetExtractor2InfoAK(file1,file2):
    '''
    We use this to extract the combinations that break AK for a particular matroid for use in search for even better bounds.

    file1, file2 (file) two .txt files that contain the sets that gave bounds in a previous run. Both files are the same.
    '''
    emlist = list()
    for line in file1:
        line = line.strip()
        se = re.findall('^S.+ sets ([{][0-9, ]+[}])',line)  # First set
        de = re.findall('^Solution.+([{][0-9, ]+[}])',line)  # Last set
        mide = re.findall('^Solution sets.+, ([{][0-9, ]+[}])',line)  # Second set
        if se == [] and mide == []:
            se = re.findall('^Solution sets ([\\[][0-9, ]+[]])',line)  # First set
            mide = re.findall('^Solution sets.+, ([\\[][0-9, ]+[]])',line)  # Second set
        if len(se) > 0 and len(de) > 0 and len(mide) > 0:
            selist = list()
            delist = list()
            midelist = list()
            for i in range(1,len(se[0]),3):
                fe = int(se[0][i])
                selist.append(fe)
            for i in range(1,len(de[0]),3):
                fe1 = int(de[0][i])
                delist.append(fe1)
            for i in range(1,len(mide[0]),3):
                fe11 = int(mide[0][i])
                midelist.append(fe11)
            emlist.append(selist)
            emlist.append(delist)
            emlist.append(midelist)
    
    if emlist == []:
        for line in file2:
            line = line.strip()
            se = re.findall('^Solution sets set\(([\\[][0-9, ]+[]])',line)  # First set
            de = re.findall('^Solution.+([\\[][0-9, ]+[]])',line)  # Last set
            mide = re.findall('^Solution sets set.+, set\(([\\[][0-9, ]+[]])',line)  # Second set
            #Solution sets set([0, 3]), set([2, 5]) and set([1, 6]) have optimal value 1.125.  # For AK
            if se == [] and mide == []:
                se = re.findall('^Solution sets ([\\[][0-9, ]+[]])',line)  # First set
                mide = re.findall('^Solution sets.+, ([\\[][0-9, ]+[]])',line)  # Second set
            if len(se) > 0 and len(de) > 0 and len(mide) > 0:
                selist = list()
                delist = list()
                midelist = list()
                for i in range(1,len(se[0]),3):
                    fe = int(se[0][i])
                    selist.append(fe)
                for i in range(1,len(de[0]),3):
                    fe1 = int(de[0][i])
                    delist.append(fe1)
                for i in range(1,len(mide[0]),3):
                    fe11 = int(mide[0][i])
                    midelist.append(fe11)
                emlist.append(selist)
                emlist.append(delist)
                emlist.append(midelist)
    return emlist