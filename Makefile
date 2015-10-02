# Short Makefile to compile "root macros"

CXX      = g++ -std=c++11
LINKER   = g++ -o
TARGET = run 

SRCDIR  = src
INCDIR  = include
OBJDIR  = obj
BINDIR  = bin
PYDIR = python

SOURCES  := $(wildcard $(SRCDIR)/*.cxx) 
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(patsubst $(SRCDIR)%.cxx,$(OBJDIR)%.o,$(SOURCES))

#get staff you need from your cmssw enviroment
scramtag = $(shell cd $$CMSSW_BASE; scram tool tag $(1) $(2))

ROOTLIBS :=  $(shell root-config --evelibs)  
ROOTFLAGS := $(shell root-config --cflags)
ROOTINC :=  -I$(call scramtag, root_interface, INCLUDE)

ROOFITLIBS :=  -L$(call scramtag, roofitcore, LIBDIR) -lRooFit -lRooFitCore -lMinuit 
ROOFITFLAGS := 
ROOFITINC :=  -I$(call scramtag, roofitcore, INCLUDE)  

PYTHON_CFLAGS := $(shell python2.6-config --cflags) 
PYTHONLIBS := $(shell python2.6-config --ldflags)
PYTHON_INC := -I$(shell python2.6-config --includes)
PYTHON_LIBPATH := -L$(shell python2.6-config --prefix)/lib

BOOST_INC := -I$(call scramtag,boost,INCLUDE)
BOOST_LIB := $(call scramtag,boost,LIB) $(call scramtag,boost_filesystem,LIB) $(call scramtag,boost_regex,LIB) $(call scramtag,boost_python,LIB) boost_iostreams
BOOST_LIB := $(patsubst %,-l%,$(BOOST_LIB)) -L$(call scramtag,boost,LIBDIR)

USERLDFLAGS += $(BOOST_LIB) $(PYTHON_LIBPATH) $(PYTHONLIBS) $(ROOTLIBS) $(ROOFITLIBS)  
USERCXXFLAGS += -g -Wall -O2 
USERCXXFLAGS += $(PYTHON_CFLAGS) $(ROOTFLAGS) $(ROOFITFLAGS)
USERINCFLAGS += $(BOOST_INC) $(PYTHON_INC) $(ROOTINC) $(ROOFITINC) -I$(INCDIR)


$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)	
	$(LINKER) $@ $(USERLDFLAGS) $(OBJECTS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	@mkdir -p $(OBJDIR)	
	$(CXX) $< $(USERINCFLAGS) $(USERCXXFLAGS) -c -o $(OBJDIR)/$(notdir $@) 
	@echo "Compiled "$<" successfully!"


.PHONEY: clean
clean:
	@rm $(OBJDIR)/*
	@echo "Cleanup complete!"

.PHONEY: distclean
distclean: clean
	@rm $(BINDIR)/*
	@echo "Executable removed!"

