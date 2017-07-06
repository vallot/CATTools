import os, time, socket, sys

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


print "Directory: " + str(sys.argv[1])
Drt  = str(sys.argv[1])

LogDrt = Drt + "_" + time.strftime('%Hh%Mm%Ss') + ".info"
os.system("find " + Drt + " -name '*.root' -size -1k -exec ls {} \+ > " + LogDrt)

fr = open(LogDrt,'r')

NSamples = 0
FileName = []
for line in fr:
    if ".root" in line:
        NSamples += 1
        sbeg = Drt + "/job_"
        send = "_cfg.py" 
        snum = sbeg + find_between( line, "Tree_ttbbLepJets_",".root" ) + send
        num = int(find_between( snum, sbeg,send ))
        if num < 100 and num > 9:
            snum = Drt + "/job_0" + str(num) + send 
        elif num < 10:
            snum = Drt + "/job_00" + str(num) + send 
        FileName.append(str(snum))
print "Number of sample to be processed -> " + str(NSamples)
fr.close()
os.system("rm -rf " + LogDrt)

OutList = open('MissingFiles_' + Drt + '.txt', 'w')
for nfile in range(0,NSamples):
    fc_name = FileName[nfile]
    if os.path.exists(fc_name):
        fc = open(fc_name,'r')
        for line in fc:
            if "sdfarm" in line:
                            #print str(line)
                OutList.write(find_between( line, "'", "'" ))
                OutList.write('\n')
        fc.close()
    else:
        print "File " + str(fc_name) + " not found...."
OutList.close()
