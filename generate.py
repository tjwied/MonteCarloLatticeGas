"""
Generate WHAM input file 
"""

import re
 
window = open('windows', 'r').readlines()
x = window[0]
print x[0] + '_gas_1d_wham.pmf'
print x[0] + '_gas_1d_wham.rho'
print x[0] + '_gas_1d_wham.bia'
print x[0] + '_gas_1d_wham.fff'
print """0 1 0.01 0
101  50000 100
0.000001"""
window = open('windows', 'r').readlines()
for file in window:
	matchobj = re.match('(\d+)\.(\d+)_k100.dat', file)
	Target1 = matchobj.group(1)
	Target2 = matchobj.group(2)
	Target3  = Target1+"."+Target2
	Target4 = float(Target3)
	k = 
	print file.strip()[0:]
	print Target3 + "    " + str(k)
