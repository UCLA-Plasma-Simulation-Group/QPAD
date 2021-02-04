# Please make sure the input file does not have any comments when using this script 

import json
from collections import OrderedDict
import getopt
import sys
import helper

# check input arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], 'ha:', ['help', 'angle='])
except getopt.GetoptError:
    print('Usage: merge_mode -a <angle in degree>')
    sys.exit(2)

# read input arguments
phi_degree = 0
for opt, arg in opts:
    if opt in ('-h', '--help'):
        print('Usage: merge_mode -a <angle in degree>')
    elif opt in ('-a', '--angle'):
        phi_degree = arg

# Find all the quantities that need to merge modes from the input file.
# Then collect their paths in the output.

with open('../qpinput.json') as finput:
    qpinput = json.load(finput,object_pairs_hook=OrderedDict)

all_paths = []

# get all the paths for the fields
fields = qpinput["field"]["diag"][0]["name"]
for field in fields:
    prefix = field.split('_')[0]
    prefix = prefix.capitalize()
    all_paths.append('../Fields/' + prefix + '/')

# get all the paths for the species
nspecies = qpinput["simulation"]["nspecies"]
if nspecies > 0:
    for i in range(nspecies):
        diags = qpinput["species"][i]["diag"]
        for j in range(len(diags)):
            names = diags[j]["name"]
            if "charge_cyl_m" in names:
                all_paths.append('../Species' + str(i+1) + '/Charge/')

# get all the paths for the neutrals
nneutrals = qpinput["simulation"]["nneutrals"]
if nneutrals > 0:
    for i in range(nneutrals):
        diags = qpinput["neutrals"][i]["diag"]
        for j in range(len(diags)):
            names = diags[j]["name"]
            if "charge_cyl_m" in names:
                all_paths.append('../Neutral' + str(i+1) + '/Charge/')
            if "ion_cyl_m" in names:
                all_paths.append('../Neutral' + str(i+1) + '/Ion_charge/')

# get all the paths for the beams
nbeams = qpinput["simulation"]["nbeams"]
if nbeams > 0:
    for i in range(nbeams):
        diags = qpinput["beam"][i]["diag"]
        for j in range(len(diags)):
            names = diags[j]["name"]
            if "charge_cyl_m" in names:
                all_paths.append('../Beam' + str(i+1) + '/Charge/')

# Merge modes in all the paths 
for path in all_paths:
    helper.merge(path,phi_degree)
