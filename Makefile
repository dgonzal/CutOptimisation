
ROOTLIBS := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags)

PYTHON_CFLAGS := $(shell python2.6-config --cflags) 
PYTHONLIBS := $(shell python2.6-config --ldflags)
PYTHON_INC := -I$(shell python2.6-config --includes)
PYTHON_LIBPATH := -L$(shell python2.6-config --prefix)/lib


USERLDFLAGS += $(PYTHON_LIBPATH) $(PYTHONLIBS) $(ROOTLIBS) 
USERCXXFLAGS := -g -Wall -O2 
USERCXXFLAGS += $(PYTHON_CFLAGS) $(ROOTFLAGS)
#USERCXXFLAGS := 

scramtag = $(shell cd $$CMSSW_BASE; scram tool tag $(1) $(2))

BOOST_INC := -I$(call scramtag,boost,INCLUDE)
BOOST_LIB := $(call scramtag,boost,LIB) $(call scramtag,boost_filesystem,LIB) $(call scramtag,boost_regex,LIB) $(call scramtag,boost_python,LIB) boost_iostreams
BOOST_LIB := $(patsubst %,-l%,$(BOOST_LIB)) -L$(call scramtag,boost,LIBDIR)


main: cutopt.o 
	g++ cutopt.o -o cutopt $(USERLDFLAGS) $(BOOST_LIB)  


cutopt.o: cutopt.cxx root-utils.hpp 
	g++  $(USERCXXFLAGS) -c cutopt.cxx $(BOOST_INC) $(PYTHON_INC)




.PHONY: clean
clean:
	rm -rf *o cutopt

