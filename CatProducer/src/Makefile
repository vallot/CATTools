ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit

# Linux with egcs
DEFINES 	= -DNO_ORCA_CLASSES
CXX		= g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD		= g++
LDFLAGS 	= -g
SOFLAGS 	= -shared

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= $(ROOTLIBS) -lEG
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES 	= $(wildcard Cat*.cc)
HEADERS 	= $(wildcard ../interface/Cat*.h)
OBJECTS 	= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS 	= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))


#PROGRAMS	=  $(SOBJECTS)

all: libToto.so clean ; cp libToto.so ../test; cp libToto.so ../skimmer

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(DEPENDS) *Dict.* core

.SUFFIXES: .$(SrcSuf) .C .o .so

###

Dict.$(SrcSuf): $(HEADERS) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@$(ROOTSYS)/bin/rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERS) ./LinkDef.h

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

libToto.so: $(OBJECTS) Dict.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@
	@echo "$@ done"

