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
DEFINES       = -DNO_ORCA_CLASSES
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD		= g++
LDFLAGS		= -g -O -Wall -fPIC
SOFLAGS		= -shared

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= $(ROOTLIBS)  -lEG 
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard ./tinyxml/*.cc)
HEADERS		= $(wildcard ./tinyxml/*.h)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

all:  libTinyXML.so;

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(DEPENDS) *Dict.* core 

.SUFFIXES: .$(SrcSuf) .C .o .so

libTinyXML.so: $(OBJECTS) 
	@echo "Building libTinyXML..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@
