#!/usr/bin/env python

for relPt in range(50,1000,50) :
  for dR in range(0,500,25) :

    print """arguments  = %f %f
queue 1 """%(relPt/1000.,dR/1000.)
