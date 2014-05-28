EMBasins.mexmaci64 : EMBasins.cpp EMBasins.h BasinModel.h TreeBasin.h BasinModel.o TreeBasin.o
	/Applications/MATLAB_R2012a.app/bin/mex -largeArrayDims  -I/usr/local/include/boost_1_52_0 -lgsl -lc++ EMBasins.cpp BasinModel.o TreeBasin.o

TreeBasin.o : TreeBasin.cpp TreeBasin.h BasinModel.h EMBasins.h
	g++ -O3 -c -I/usr/local/include/boost_1_52_0 -std=c++11 -stdlib=libc++ TreeBasin.cpp

BasinModel.o : BasinModel.cpp BasinModel.h EMBasins.h
	g++ -O3 -c -stdlib=libc++ -std=c++11 BasinModel.cpp
