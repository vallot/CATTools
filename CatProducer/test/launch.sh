#!/bin/bash

usage(){
  echo "####################################################################################################"
  echo "              LAUNCH GRID JOB"
  echo ""
  echo "The first test consist to write a our jdl file (see example: template.jdl) or use makeJDL command"
  echo "Second step: having a proxy active"
  echo ""
  echo "usage::"
  echo "-makeJDL 			myCode myRootOutputFile.root  -> will create a test.jdl"
  echo "-test myJDLfile.jdl: 		check if jdl file is conform"
  echo "-submit myJDLfile.jdl myLabel:  will use glite-wms2.conf file to submit the job"
  echo "-status myLabel:		... "
  echo "-getoutput myLabel:		which retrieve the output and store it in $PWD/myLabel"
  echo ""
  echo "Warning:: myLabel should be the same variable for all the steps"
  echo "####################################################################################################"
}

if [ $# -lt 2 ] ; then
   echo "missing argument"
   usage
fi
				  

  case "$1" in
  -makeJDL)
          export var=`pwd`
          echo $var
	  sed -e 's/EXECUTABLE/'$2'/g' template.jdl | sed -e 's/OUTPUT/'$3'/g' > test.jdl
	  sed -e 's/EXECUTABLE/'$2'/g' template.sh  > test.sh
	  ;;
  -test)   
         if [ $# -eq 2 ]; then
           if [ -e  $2 ]; then
             if [ -e glite-wms2.conf ]; then
	       echo glite-wms-job-list-match -a -c $2
	       glite-wms-job-list-match -a -c glite-wms2.conf $2
              else
	        echo "glite-wms2.conf doesn't exist"
	      fi
	   else
	     echo $2 "doesn't exist"
	   fi
	 else
	   echo  "missing argument"
	   usage
	 fi
	 ;;
  -submit) 
         if [ $# -eq 3 ]; then
           if [ -e  $2 ]; then
             if [ -e glite-wms2.conf ]; then
                echo glite-wms-job-submit -a -c glite-wms2.conf -o ./.$3 $2
                glite-wms-job-submit -a -c glite-wms2.conf -o ./.$3 $2
	     else 
	        echo "glite-wms2.conf doesn't exist"
	     fi
	   else
	     echo $2 "doesn't exist"
           fi
	 else
	   echo  "missing argument"
	   usage
	 fi
	 ;;
  -status) 
         if [ $# -eq 2 ]; then
            echo glite-wms-job-status -i .$2
            glite-wms-job-status .$2
	 else
	   echo  "missing argument"
	   usage
	 fi
	 ;;
  -getoutput) 
         if [ $# -eq 2 ]; then
            echo glite-wms-job-output --dir $PWD/$2 -i .$2
            glite-wms-job-output --dir $PWD/$2 -i .$2
	 else
	   echo  "missing argument"
	   usage
	 fi
	 ;;
  *)
    usage
    ;;					      
  esac
 shift
	 
		 
