# Global variables 
from gurobipy import *

## For an 8-point matroid, we would have:
## Part = 8 (one member represents the dealer while the other 7 are the participants)
## vrbls = 9 (one auxiliary variable is added)
## The number of auxiliary variables added depends on the depth we're checking.
## If for 1-CI/AK, we add just one, while 3 would be added for 3-CI/AK.

vrbls = 10  # participants + dealer + auxiliary variable(s) 
Part = 9  # participants + dealer

p = Model("gurotest")
w = p.addVars(range(0,2**vrbls+1),name="w")


destiny = {
' Tictactoe port 5 '  :  ['12346', '12347', '12367', '12368', '12378', '12467', '12468', '12478', '13467', '13468', '13478', '13678', '23467', '2348', '2678', '34678', '51234', '51237', '51246', '51267', '5128', '51347', '51348', '51367', '51368', '51378', '51467', '51468', '51478', '51678', '52347', '52367', '52368', '52378', '52467', '52468', '52478', '5346', '53478', '54678']
}
'''
destiny = {
' Bas_F8 port 1 '  :  ['134', '245', '127', '236', '156', '1467', '1357', '2467', '2357', '3467', '3457', '4567', '3567'] ,
' Bas_V8 port 3 '  :  ['312', '245', '3246', '3247', '3256', '3267', '3257', '3146', '3147', '3156', '3157', '3456', '3467', '3457', '3567', '1246', '1247', '1256', '1267', '1257', '2467', '2567', '1456', '1467', '1457', '1567'] ,
' Bas_V8 port 1 '  :  ['123', '145', '1246', '1247', '1256', '167', '1257', '1346', '1347', '1356', '1357', '2346', '2347', '2456', '2467', '2457', '2356', '2367', '2357', '2567', '3456', '3467', '3457', '3567'] 
}
'''