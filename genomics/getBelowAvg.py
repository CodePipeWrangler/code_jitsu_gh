#!/usr/bin/env python
# NAME
#   getBelowAvg.py - Filter *.coords files to hits with below avg % IDY.
#
# SYNOPSIS
#   ./top_line.awk []
#   runs in theworking directory
#
# OPERANDS
#     INPUT_FILE
#
#
# EQUIVALENT TO
#

# Method to list files in directory containing {""}.format(*coords)...
import glob # finds all the pathnames matching a specified pattern
allCoords = list(glob.glob('*coords'))

import pandas as pd # data analysis and manipulation tools

#stat = dict.fromkeys(['min', 'max', 'avg']) #not needed but creates empty dictionary with keys

# The header originally does not have enough column names. I created 'TAGS1(2)' from 'TAGS'.
header = ['[S1]', '[E1]', '[S2]', '[E2]', '[LEN 1]', '[LEN 2]', '[% IDY]',
       '[LEN R]', '[LEN Q]', '[COV R]', '[COV Q]', '[TAGS1]','[TAGS2]']

# get stats and rows containing below average % IDY for all coords files in directory.
for file in allCoords:
        new = pd.read_csv(file,sep='\t',skiprows=(0,1,2)) # load file
        # tidy up df
        new.reset_index(inplace=True) # move 1st column from index
        new.columns=header # place in header
        # get stats
        stat = {'max':new.iloc[:,6].max(),'min':new.iloc[:,6].min(),'avg':new.iloc[:,6].mean()}
        print(file,'\t',stat) #print stats to screen
        worthy = new.loc[new['[% IDY]'] <= stat['avg']] # get rows containing % IDY below avg
        worthy.to_csv('IDYv{:0.3f}_{}'.format(stat['avg'],file), sep='\t') # save as tab-delimitted file
                # note '0.Xf' sets the decimal places. The {} are used to insert values
del stat # remove the statistical variables
