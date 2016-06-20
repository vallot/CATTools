import os,commands
import sys

if len(sys.argv) is 1:
  print "Please, add the name as like followings."
  print "> python job_check.py .  \n"
  sys.exit()

arg = sys.argv[1]
xrootdpath = "/xrootd/store/user//youngjo/Cattools/v7-6-5v1/"
v = "v1" 



##########################################
##########################################
##########################################
def lsSub(path):
  a= len([f for f in os.listdir(path) if f.find("cfg")>-1 ])
  b= len([f for f in os.listdir(path) if f.find("err")>-1 ])
  return a,b

def ls(path):
  a=[f for f in os.listdir(path) if  os.path.isdir(os.path.join(path,f)) ]
  cc = {}
  for b in a:
    aaa,bbb = lsSub(path+"/"+b)
    #if bbb is 0 :
    cc[b] = aaa,bbb

  return cc

def lsXrootdSub(path):
  a= len([f for f in os.listdir(path) if f.find(".root")>-1 ])
  return a

def lsXrootd(path,v,ccc):
  a=[f for f in os.listdir(path) if  os.path.isdir(os.path.join(path,f)) and f.find(v)>-1 ]
  cc = {}
  for b in a:
    aaa = lsXrootdSub(path+"/"+b)
    bb = b.replace(v,"")
    if not bb in ccc.keys():
      if aaa is 0 : 
        cc[bb] = 0,0,aaa 
    elif ccc[bb][0] != aaa or ccc[bb][0] != ccc[bb][1] :
      cc[bb] = ccc[bb][0],ccc[bb][1],aaa

  return cc



##########################################
##########################################
##########################################
##########################################
dd=ls(arg)
#print "  "+str(dd)

ee=lsXrootd(xrootdpath,v,dd)
print "  "+str(ee)

from sets import Set
ff = Set()
for gg in ee.keys():
  ff.add(gg)

print "ff : "+str(ff)



##########################################
##########################################
##########################################


