
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

cp ../src/libToto.so .

g++  -L `pwd` -l Toto -I `root-config --incdir` `root-config --libs` catContentDump.C -o catContentDump
	
if [ $? -eq 0 ]
then
	echo "Compiling OK, now run the code!"

	if [ ! $# -eq 2 ]; then
	  exit;
 
        fi

        file=""
	
	if [ $1 == "dcap" ];then

	  for i in $(ls $2 | grep -v BADFILE);do

	    file="$file;dcap://maite.iihe.ac.be$2/$i"

	  done

        fi

	if [ $1 == "local" ];then

	  for i in $(ls $2);do

	    file="$file;$2/$i"

	  done

        fi

	#echo $file

	./catContentDump --inputfiles "$file"
else
	echo "Problem while compiling!"
fi
