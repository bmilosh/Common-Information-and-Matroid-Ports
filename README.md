# Common-Information-and-Matroid-Ports
These are files corresponding to results of the paper "Common Information, Matroid Representation, and Secret Sharing for Matroid Ports",
joint work by Michael Bamiloshin, Aner Ben-Efraim, Oriol Farras, and Carles Padro.

They contain codes to check if a matroid is 1-CI or 1-AK, using LP techniques.
They also contain codes to check for bounds on lambda and sigma for a matroid.

AKextensions.py is used to check if a matroid is 1-AK. Can be updated to specify the kind of subsets to check.
CIextensions.py is used to check if a matroid is 1-CI. Can be updated to specify the kind of subsets to check.

checkSigma.py is used to find bounds on sigma for a given matroid. Can be updated to specify the kind of subsets to check.
checkSigma_Improved.py is used to improve bounds on sigma derived using checkSigma.py

checkLambda.py is used to find bounds on lambda for a given matroid. Can be updated to specify the kind of subsets to check.
checkLambda_Improved.py is used to improve bounds on lambda derived using checkLambda.py

fibonew2.py contains all functions needed to run CIextensions.py, checkSigma.py, e.t.c. 

config.py hosts the global variables Part and vrbls. These MUST be updated whenever the parameters change 
(size of the matroid's ground set and/or the number of auxiliary variables used).

circuitgenerator.py contains functions to generate the matroid ports. 
Since some ports of a particular matroid might be isomorphic, this can be checked offline using SageMath (or whatever resource the user prefers).

mat84bases.py contains ALL rank-4 matroids on 8 points. 
The resource "Matroids on 9 Elements" provided by Gordon Royle and Dillon Mayhew 
(https://research-repository.uwa.edu.au/en/datasets/matroids-on-9-elements) contains ALL matroids on up to 9 points.
The matroid ids in our paper correspond to their numbers as they appear in this Royle and Mayhew database, unless otherwise stated.

All files were prepared with python 3.7.4, while our LP solver is Gurobi version 9.0.0

test_file.ipnyb is a Jupyter notebook that contains some example runs.
