##
## Velour -- memory efficient short read de novo DNA sequence assembly tool for
##           gigabase-scale genomes
##
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#
##  MAKE TARGET  |   BINARY     |  DESCRIPTION
## --------------|--------------|----------------------------------------------
##  bench        | velour.bench | for timing runs: fully optimized, no asserts, no debug info
## opt (default) | velour       | for general use: optimized, most asserts, some debug info
##  debug        | velour.debug | for debugger use: no optimization, all asserts, verbose debug info
##
##
## Building strategies:
##
## * build in the base directory:
##     localhost$ make
##
## * build in an arbitrary directory:
##     localhost$ cd builddir && make -C basedir
##

MAXKMERLENGTH=31

OSNAME=$(shell uname -s)
ifeq ($(OSNAME),Darwin)
  ifeq ($(shell sysctl -n hw.optional.x86_64),1)
    ARCH=$(shell if [ "`sysctl -n hw.memsize`" -gt 4294967296 ] ; then echo x86_64 ; else echo i386 ; fi)
  else
    ARCH=i386
  endif
else
  ARCH=$(shell uname -p)
endif

##
## partitioner selection
##   uncomment the target partitioner and its related variables
##

#PARTITIONER=mk_random.cpp
#PARTITIONER_MIDDLE_MINIKMER='true'

#PARTITIONER=mk_scooping.cpp
#PARTITIONER_MIDDLE_MINIKMER='true'

PARTITIONER=mk_static.cpp
PARTITIONER_POWEROF4='true'
PARTITIONER_MIDDLE_MINIKMER='true'

#PARTITIONER=mk_greedy_clustering.cpp
#PARTITIONER_TRAINS='true'
#PARTITIONER_MINIKMER_BITS='true'

##
## END partitioner selection
##


## default compiler
##   to override: make CXX=g++-4.2
ifeq ($(OSNAME),Darwin)
  CXX=$(shell if [ ! -z `which g++-4.2` ] ; then echo g++-4.2 ; else echo g++ ; fi)
  #CXX=icpc
else
  CXX=g++
  #CXX=icpc
endif

## generic compilation flags
##     to override: make CXXFLAGS="-Wall"
##   to supplement: make ADD_CXXFLAGS="-Wnone" ADD_LDFLAGS="-lm"
##
COMMON_CXXFLAGS =
ifeq ($(ARCH),x86_64)
  COMMON_CXXFLAGS += -m64 -DARCH_64BIT
else
  COMMON_CXXFLAGS += -m32 -DARCH_32BIT
endif

BENCH_CXXFLAGS = $(COMMON_CXXFLAGS) -DNDEBUG
OPT_CXXFLAGS   = $(COMMON_CXXFLAGS)
DEBUG_CXXFLAGS = $(COMMON_CXXFLAGS)

COMMON_LDFLAGS = -lm
ifeq ($(ARCH),x86_64)
  COMMON_LDFLAGS += -m64
else
  COMMON_LDFLAGS += -m32
endif

BENCH_LDFLAGS = $(COMMON_LDFLAGS)
OPT_LDFLAGS   = $(COMMON_LDFLAGS)
DEBUG_LDFLAGS = $(COMMON_LDFLAGS)

BENCH_CXXFLAGS += $(ADD_CXXFLAGS)
OPT_CXXFLAGS   += $(ADD_CXXFLAGS)
DEBUG_CXXFLAGS += $(ADD_CXXFLAGS)

BENCH_LDFLAGS += $(ADD_LDFLAGS)
OPT_LDFLAGS   += $(ADD_LDFLAGS)
DEBUG_LDFLAGS += $(ADD_LDFLAGS)

## g++ specific flags
##
ifneq ($(findstring g++,$(CXX)),)
  RECOGNIZED = "g++"

  S_CXXFLAGS = -Wall -march=native -fno-strict-aliasing
  S_CXX_BENCH = $(S_CXXFLAGS) -O3
  S_CXX_OPT   = $(S_CXXFLAGS) -O3 -g
  S_CXX_DEBUG = $(S_CXXFLAGS) -O0 -g3 -ggdb3 -fno-omit-frame-pointer

  S_LDFLAGS =
  S_LD_BENCH = $(S_LDFLAGS)
  S_LD_OPT   = $(S_LDFLAGS)
  S_LD_DEBUG = $(S_LDFLAGS)
endif

## Intel icc specific flags
##
##  TODO: flags to later consider: align, no-prev-div
ifneq ($(findstring icpc,$(CXX)),)
  RECOGNIZED = "icpc"

  S_CXXFLAGS = -Wall -xHost -fno-exceptions -fno-rtti
  S_CXXFLAGS += -shared-intel -shared-libgcc
  S_CXX_BENCH = $(S_CXXFLAGS) -O3 -ipo
  S_CXX_OPT   = $(S_CXXFLAGS) -O3 -g -ipo
  S_CXX_DEBUG = $(S_CXXFLAGS) -O0 -g -debug full -fno-omit-frame-pointer

  S_LDFLAGS =
  S_LD_BENCH = $(S_LDFLAGS)
  S_LD_OPT   = $(S_LDFLAGS)
  S_LD_DEBUG = $(S_LDFLAGS)
endif

ifeq ($(RECOGNIZED),)
  $(error Makefile did not recognize compiler CXX=$(CXX))
endif


## code path override
##     to override: make CODEPATH="-DFOO"
##   to supplement: make ADD_CODEPATH="-DBAR"
ifndef CODEPATH
CODEPATH += -DMAXKMERLENGTH=$(MAXKMERLENGTH)
CODEPATH += -DSMALL_NODES
#CODEPATH += -DUNIQUE
#CODEPATH += -DCOLOR_8BIT
#CODEPATH += -DVELOUR_TBB
#CODEPATH += -DVELOUR_TBB -DTBB_USE_DEBUG
#CODEPATH += -DUSE_TBB_ALLOC
#CODEPATH += -DUSE_LIBC_ALLOC
#CODEPATH += -DVERIFY
#CODEPATH += -DVERIFY_THREADSAFE
#CODEPATH += -DITT_THREAD_CHECKER
#CODEPATH += -DITT_THREAD_PROFILER
#CODEPATH += -DVALGRIND
CODEPATH += -DVELVET_EMULATION
endif
CODEPATH += $(ADD_CODEPATH)

## input partitioning
ifdef PARTITIONER
  CODEPATH += -DPARTITIONING
  ifdef PARTITIONER_TRAINS
	CODEPATH += -DPART_TRAINS_ON_INPUT
  endif
  ifdef PARTITIONER_MINIKMER_BITS
	CODEPATH += -DPART_MINIKMER_BITS
  endif
  ifdef PARTITIONER_POWEROF4
	CODEPATH += -DPART_POWEROF4
  endif
  ifdef PARTITIONER_MIDDLE_MINIKMER
	CODEPATH += -DPART_MIDDLE_MINIKMER
  endif
endif

ifneq ($(CODEPATH),)
  DEFINES += -DOVERRIDE_CODEPATH $(CODEPATH)
endif

ifneq ($(findstring TBB,$(CODEPATH)),)
BENCH_LDFLAGS += -ltbb -ltbbmalloc_proxy -ltbbmalloc -lrt
OPT_LDFLAGS   += -ltbb -ltbbmalloc_proxy -ltbbmalloc -lrt
DEBUG_LDFLAGS += -ltbb_debug -ltbbmalloc_proxy_debug -ltbbmalloc_debug -lrt
endif

ifneq ($(findstring ITT,$(CODEPATH)),)
BENCH_LDFLAGS += -littnotify
OPT_LDFLAGS   += -littnotify
DEBUG_LDFLAGS += -littnotify
endif

##
## pre-processor parameter defines
##     to override: make USERDEFINES="-DFOO -D'BAR=37'"
##   to supplement: make ADD_USERDEFINES="-DFOOBAR"
ifndef USERDEFINES
#USERDEFINES += -DVERBOSE
USERDEFINES += -D'RANDOMSEED=37'
endif
USERDEFINES += $(ADD_USERDEFINES)
DEFINES += $(USERDEFINES)



#############################################################################
#############################################################################
###
### No user-configurable variables beyond this point?
###
#############################################################################
#############################################################################

BASEDIR=$(CURDIR)
SRCDIR=$(BASEDIR)/src
OBJBASEDIR=$(PWD)/obj
BINDIR=$(PWD)
export BASEDIR
export BINDIR

VPATH=$(SRCDIR)

BENCH_OBJDIR=$(OBJBASEDIR)/bench
DEBUG_OBJDIR=$(OBJBASEDIR)/debug
OPT_OBJDIR=$(OBJBASEDIR)/opt

# list of binaries
BASE_BINARIES=velour

BENCH_BINARIES=$(addsuffix .bench, $(addprefix $(BINDIR)/,$(BASE_BINARIES)))
DEBUG_BINARIES=$(addsuffix .debug, $(addprefix $(BINDIR)/,$(BASE_BINARIES)))
OPT_BINARIES=$(addprefix $(BINDIR)/, $(BASE_BINARIES))

export VELOUR=$(BINDIR)/velour
export VELOURDEBUG=$(BINDIR)/velour.debug
export VELOURBENCH=$(BINDIR)/velour.bench

COMMON_SRCS = allocators.cpp \
        distance.cpp \
        flowing.cpp \
	globals.cpp \
        graphviz.cpp \
        kmer.cpp \
        kmerGraph.cpp \
        kmerNode.cpp \
        minikmer.cpp \
        node_allocators.cpp \
        parsing.cpp \
        partition.cpp \
        $(PARTITIONER) \
        pg_tipClipping.cpp \
        preGraphConstruction.cpp \
	pregraph_distribution.cpp \
        pregraph_partitioning.cpp \
        quilt.cpp \
        seqGraph.cpp \
        seqGraph2.cpp \
        seqNode.cpp \
        sequence.cpp \
        slicing2.cpp \
	sg_bubblePopping.cpp \
	sg_covcutoff.cpp \
        sg_tipClipping.cpp \
        sg_verify.cpp \
        split.cpp \
        stat_components.cpp \
        utility.cpp

MAIN_SRCS = velour.cpp

COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

COMMON_BENCHOBJS=$(addprefix $(BENCH_OBJDIR)/,$(COMMON_OBJS))
COMMON_DEBUGOBJS=$(addprefix $(DEBUG_OBJDIR)/,$(COMMON_OBJS))
COMMON_OPTOBJS=$(addprefix $(OPT_OBJDIR)/,$(COMMON_OBJS))

MAIN_BENCHOBJS=$(addprefix $(BENCH_OBJDIR)/,$(MAIN_OBJS))
MAIN_DEBUGOBJS=$(addprefix $(DEBUG_OBJDIR)/,$(MAIN_OBJS))
MAIN_OPTOBJS=$(addprefix $(OPT_OBJDIR)/,$(MAIN_OBJS))

.PRECIOUS: Makefile
.PHONY: default all bench debug opt clean veryclean print-gitbranch

default: opt

all: bench debug opt

bench: $(BENCH_BINARIES)
	@if [ $(MAXKMERLENGTH) -gt 31 ] ; then \
		echo -n "NOTE: MAXKMERLENGTH = $(MAXKMERLENGTH) > 31, using large-kmer data structure." ; \
		echo "  Velour run time time and memory use will be impacted." ; \
	 fi

debug: $(DEBUG_BINARIES)
	@if [ $(MAXKMERLENGTH) -gt 31 ] ; then \
		echo -n "NOTE: MAXKMERLENGTH = $(MAXKMERLENGTH) > 31, using large-kmer data structure." ; \
		echo "  Velour run time time and memory use will be impacted." ; \
	 fi

opt: $(OPT_BINARIES)
	@if [ $(MAXKMERLENGTH) -gt 31 ] ; then \
		echo -n "NOTE: MAXKMERLENGTH = $(MAXKMERLENGTH) > 31, using large-kmer data structure." ; \
		echo "  Velour run time time and memory use will be impacted." ; \
	 fi

$(notdir $(BENCH_BINARIES)):
	$(MAKE) $(BINDIR)/$@

$(notdir $(DEBUG_BINARIES)):
	$(MAKE) $(BINDIR)/$@

$(notdir $(OPT_BINARIES)):
	$(MAKE) $(BINDIR)/$@


$(filter $(BINDIR)/velour%, $(BENCH_BINARIES)): MAIN_OBJ = $(BENCH_OBJDIR)/velour.o
$(filter $(BINDIR)/velour%, $(BENCH_BINARIES)): $(BENCH_OBJDIR)/velour.o
$(filter $(BINDIR)/velour%, $(DEBUG_BINARIES)): MAIN_OBJ = $(DEBUG_OBJDIR)/velour.o
$(filter $(BINDIR)/velour%, $(DEBUG_BINARIES)): $(DEBUG_OBJDIR)/velour.o
$(filter $(BINDIR)/velour%, $(OPT_BINARIES)): MAIN_OBJ = $(OPT_OBJDIR)/velour.o
$(filter $(BINDIR)/velour%, $(OPT_BINARIES)): $(OPT_OBJDIR)/velour.o


$(BENCH_BINARIES): CXXFLAGS = $(BENCH_CXXFLAGS) $(S_CXX_BENCH) $(DEFINES) -I$(SRCDIR)
$(BENCH_BINARIES): LDFLAGS = $(BENCH_LDFLAGS) $(S_LD_BENCH) $(DEFINES)
$(BENCH_BINARIES): $(COMMON_BENCHOBJS)
	$(CXX) $(MAIN_OBJ) $(COMMON_BENCHOBJS) $(LDFLAGS) -o $@

$(DEBUG_BINARIES): CXXFLAGS = $(DEBUG_CXXFLAGS) $(S_CXX_DEBUG) $(DEFINES) -I$(SRCDIR)
$(DEBUG_BINARIES): LDFLAGS = $(DEBUG_LDFLAGS) $(S_LD_DEBUG) $(DEFINES)
$(DEBUG_BINARIES): $(COMMON_DEBUGOBJS)
	$(CXX) $(MAIN_OBJ) $(COMMON_DEBUGOBJS) $(LDFLAGS) -o $@

$(OPT_BINARIES): CXXFLAGS = $(OPT_CXXFLAGS) $(S_CXX_OPT) $(DEFINES) -I$(SRCDIR)
$(OPT_BINARIES): LDFLAGS = $(OPT_LDFLAGS) $(S_LD_OPT) $(DEFINES)
$(OPT_BINARIES): $(COMMON_OPTOBJS)
	$(CXX) $(MAIN_OBJ) $(COMMON_OPTOBJS) $(LDFLAGS) -o $@


%/dummy.txt:
	@mkdir -p $(@D)
	@touch $@

$(MAIN_BENCHOBJS) $(COMMON_BENCHOBJS): $(BENCH_OBJDIR)/dummy.txt Makefile
	$(CXX) $(CXXFLAGS) -c $(patsubst %.o,%.cpp,$(addprefix $(SRCDIR)/,$(notdir $@))) -o $@

$(MAIN_DEBUGOBJS) $(COMMON_DEBUGOBJS): $(DEBUG_OBJDIR)/dummy.txt Makefile
	$(CXX) $(CXXFLAGS) -c $(patsubst %.o,%.cpp,$(addprefix $(SRCDIR)/,$(notdir $@))) -o $@

$(MAIN_OPTOBJS) $(COMMON_OPTOBJS): $(OPT_OBJDIR)/dummy.txt Makefile
	$(CXX) $(CXXFLAGS) -c $(patsubst %.o,%.cpp,$(addprefix $(SRCDIR)/,$(notdir $@))) -o $@

clean:
	-rm -f $(BENCH_BINARIES) $(DEBUG_BINARIES) $(OPT_BINARIES)
	-if [ "$(OSNAME)" = "Darwin" ] ; then rm -rf $(addsuffix .dSYM, $(DEBUG_BINARIES) $(OPT_BINARIES)) ; fi
	-if [ -d $(OBJBASEDIR) ] ; then rm -rf $(OBJBASEDIR) ; fi
#	-rm -f $(subst $(PWD)/,,$(BENCH_BINARIES) $(DEBUG_BINARIES) $(OPT_BINARIES))

checkclean:

-include Makefile.check

veryclean: clean checkclean

# emit current git branch name, if any
print-gitbranch:
ifneq ($(wildcard $(BASEDIR)/.git),)
  ifneq ($(shell which git),)
	@-echo -n "git branch:  " && git symbolic-ref HEAD|sed s,refs/heads/,, || echo "???"
  endif
endif

#
# generate dependence files
#
BENCHDEPS=$(MAIN_BENCHOBJS:.o=.d) $(COMMON_BENCHOBJS:.o=.d)
DEBUGDEPS=$(MAIN_DEBUGOBJS:.o=.d) $(COMMON_DEBUGOBJS:.o=.d)
OPTDEPS=$(MAIN_OPTOBJS:.o=.d) $(COMMON_OPTOBJS:.o=.d)

$(BENCH_OBJDIR)/%.d: $(SRCDIR)/%.cpp $(BENCH_OBJDIR)/dummy.txt Makefile
	@set -e; rm -f $@; \
	  $(CXX) -MM -MT $(patsubst %.d,%.o,$@) $(CXXFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$

$(DEBUG_OBJDIR)/%.d: $(SRCDIR)/%.cpp $(DEBUG_OBJDIR)/dummy.txt Makefile
	@set -e; rm -f $@; \
	  $(CXX) -MM -MT $(patsubst %.d,%.o,$@) $(CXXFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$

$(OPT_OBJDIR)/%.d: $(SRCDIR)/%.cpp $(OPT_OBJDIR)/dummy.txt Makefile
	@set -e; rm -f $@; \
	  $(CXX) -MM -MT $(patsubst %.d,%.o,$@) $(CXXFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$

# no need to generate dependences when cleaning
ifneq "$(MAKECMDGOALS)" "clean"
  ifneq "$(MAKECMDGOALS)" "veryclean"
    ifneq "$(MAKECMDGOALS)" "checkclean"
-include $(BENCHDEPS)
-include $(DEBUGDEPS)
-include $(OPTDEPS)
    endif
  endif
endif
