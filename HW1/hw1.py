#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 18:01:03 2017

@author: mengranzhang
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:27:37 2017

@author: mengranzhang
"""
#import sys
#result = open(sys.argv[1],'r')
#datafile = open(sys.argv[2],'r')
result = open('little.log')
datafile = open('little.dat')
#old portfolio

lines = datafile.readlines();
datafile.close()

#print lines[0]
firstline = lines[0].split()
#print "first line is", firstline

numsec = int(firstline[1])
numscen = int(firstline[3])
r = float(firstline[5])
#print "\n"
#print "number of securities:", numsec,"number of scenarios", numscen,"r",r
#print "\n"

#result from Q1
newlines = result.readlines()
result.close()

position= [0] * (numsec +1)
for i in range(0,numsec +1):
    items = newlines[i].split()
    position[i] = float(items[2]) 



#allocate prices as one-dim array
p = [0]*(1 + numsec)*(1 + numscen)


# line k+1 has scenario k (0 = today)

import numpy 

allscore = [0]*100
for i in range(0,100):
    
    k = 0
# line k+1 has scenario k (0 = today)
    while k <= numscen:
        thisline = lines[k + 1].split()
        #print "line number", k+1,"is", thisline
    # should check that the line contains numsec + 1 words
        j = 1
        #print "scenario", k,":"
        p[k*(1 + numsec)] = 1 + r*(k != 0)
        while j <= numsec:
            value = float(thisline[j])
            p[k*(1 + numsec) + j] = value * (1 + numpy.random.uniform(-0.05, 0.05))
           # print " sec ", j, " -> ", p[k*(1 + numsec) + j]
            j += 1
        k += 1
    
    pm = numpy.asarray(p).reshape(numscen +1 , numsec +1)
    final = numpy.dot(pm, numpy.asarray(position))

    score = len(final[final>0])
    allscore[i] = score
    
print allscore

from matplotlib import pyplot

pyplot.hist(allscore)
