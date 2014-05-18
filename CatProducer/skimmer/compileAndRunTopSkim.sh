
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

#cp ../src/libToto.so .

g++  -L `pwd` -l TinyXML -l Toto -I `root-config --incdir` `root-config --libs` TopSkim.C -o TopSkim

if [ $? -eq 0 ]
then
	echo "Compiling OK, now run the code!"
	./TopSkim
else
	echo "Problem while compiling!"
fi
