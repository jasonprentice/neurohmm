CC = g++
LOCAL_INCLUDE = ./include
BOOST_INCLUDE = /usr/local/include/boost_1_52_0
INCLUDE = $(LOCAL_INCLUDE):$(BOOST_INCLUDE)
CPPFLAGS = -O3 -Wall -Werror -std=c++11 -stdlib=libc++ -I$(INCLUDE)
OBJS = BasinModel.o IndependentBasin.o TreeBasin.o RNG.o paramsStruct.o 
LOCAL_LIB = ./lib
HEADERS = $(LOCAL_INCLUDE)/EMBasins.h $(LOCAL_INCLUDE)/HMM.h $(LOCAL_INCLUDE)/paramsStruct.h

FitBasinModel.mexmaci64 : FitBasinModel.cpp $(HEADERS)  $(LOCAL_LIB)/libembasins.a
	@/Applications/MATLAB_R2012a.app/bin/mex -largeArrayDims -L$(LOCAL_LIB) -I$(LOCAL_INCLUDE) FitBasinModel.cpp -lembasins -lgsl
	@mv FitBasinModel.mexmaci64 ../

libembasins.a : $(OBJS)
	@ar cr libembasins.a $(OBJS)
	@mv libembasins.a $(LOCAL_LIB)

-include $(OBJS:.o=.d)

%.o : %.cpp
	@$(CC) -c $(CPPFLAGS) $<
	@$(CC) -MM $(CPPFLAGS) $*.cpp > $*.d

clean : 
	@rm -f *.o *.d
