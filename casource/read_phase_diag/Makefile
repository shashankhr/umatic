.SUFFIXES:
.SUFFIXES: .c .o .d
vpath %.c ..

%.o: %.c	
	$(CC) $(CFLAGS) -c $<

#DEPEND rule copied from gnu make manual
%.d: %.c
	set -e; cc -MM $(CFLAGS) $< \
| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
[ -s $@ ] || rm -f $@

#for all systems
REG_LCDEFS = -Wall -DSYSNAME=$(SYSNAME) -DBL_COMPRESS 
EXE = 
DBMFLAGS = -DDMALLOC
DFLAGS= -DDEBUG_TERMS

#CYGWIN (windows pc) section
ifneq (,$(findstring CYGWIN,$(ARCH)))
WOFF = 
#DFLAGS = -DALL_SIGS
SYSLCOPTS =   
DFLAGS =
CC = gcc
REG_OPT = -O3 -ffast-math 
DEB_OPT = -g
PAR = 2
LLDLIBS = -lm 
LDFLAGS = -static -s -L/usr/local/lib $(LLDLIBS)
DEB_LIB = 
SYSNAME = cygwin
NAME = cyg
EXE = .exe

endif

#LINUX (pc) section
#also defaults not changed by sections below
ifneq (,$(findstring Linux,$(ARCH)))
WOFF = 
#DFLAGS = -DALL_SIGS
DFLAGS =
CC =  gcc
LINT = splint -weak
REG_OPT =  -ffast-math #-march=pentium3
DEB_OPT = -ggdb
PAR =  -j6

#LLDLIBS =   -L/usr/local/lib -lm -lgd
LLDLIBS =  libreadlist.a -L/usr/local/lib -lm -lgd
DEB_LIB = #-lefence 
DBM_LIB = -ldmalloc

NAME = linux
SYSNAME = linux
SYSLCOPTS =  -DLINUX
endif

#
#specific for hive (with different location for gcc3)
#
ifneq (,$(findstring hive,$(HOSTNAME)))
CC =  gccv333
#CC = icc
REG_OPT = # -O0 -xW -ipo
#gnu profiler options
#should not use debugging malloc libraray with profiling!
#REG_OPT = -g -pg  
DEB_OPT = -ggdb
NAME = gcc3
endif

#
#specific for Intel ICC complier
#looks in the path for the intel files
#
#ifneq (,$(findstring intel,$(PATH)))
#CC = icc
#REG_OPT =  -O3 -xW -ipo  -parallel # -tpp7 -openmp
#DEB_OPT = -g
#NAME = icc
#endif

#####
#end of system specific sections
#
####
MSG = "Finished making"

ifeq ($(strip $(ARCH)),)
CC = echo
SYSNAME = ERROR_NO_ARCH
NAME = ERROR_NO_ARCH
MSG = "ERROR ERROR ERROR !!!!! -- no ARCH variable detected! Please set the environment variable to describe the system. Quitting with error: "
endif

SOURCES1 =\
readphasediag.c\
phasediag.c\
readwrite_compress.c\
ph_input_struct_reader.c\
writephasediag.c

SOURCES2 =\
convert.c\
freecsv.c\
init.c\
readwrite.c\
safeopen.c\

OBJ = $(SOURCES1:.c=.o)
OBJ2 = $(SOURCES2:.c=.o)

ALLSOURCES = $(SOURCES1) $(SOURCES2)

DEPEND = $(ALLSOURCES:.c=.d)
PROG = phasediag_$(NAME)$(EXE)
RPROG = readph_$(NAME)$(EXE)
RDPROG = readph_$(NAME)deb$(EXE)
DPROG = phasediag_$(NAME)deb$(EXE)

CFLAGS = $(OPTIMIZER) $(POPTS) $(LCOPTS) $(REG_LCDEFS)


default: $(PROG)
debug: clean $(DPROG) 
read: clean $(RPROG)
new: clean $(PROG)
rdebug: clean $(RDPROG)


$(PROG) $(DPROG): POPTS = -DTESTMAIN  -DCA_DOUBLE
$(PROG) $(DPROG) $(RPROG) $(RDPROG): CFILES = $(SOURCES1) 
$(RPROG) $(RDPROG):  POPTS = -DREADMAIN -DCA_DOUBLE


$(PROG) $(RPROG) : OPTIMIZER = $(REG_OPT)
$(PROG) $(RPROG) : LCOPTS = $(SYSLCOPTS)
$(PROG) $(RPROG) : LLDFLAGS = $(LLDLIBS)

$(DPROG) $(RDPROG) : OPTIMIZER = $(DEB_OPT)
$(DPROG) $(RDPROG) : LCOPTS = $(SYSLCOPTS) $(DFLAGS)
$(DPROG) $(RDPROG) : LLDFLAGS = $(LLDLIBS) $(DEB_LIB)

include $(DEPEND)




$(PROG) $(DPROG) $(RPROG) $(RDPROG) : $(OBJ2) $(OBJ)	$(DEPEND)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(OBJ2) $(LLDFLAGS);echo "Finished making" $@


clean:
	rm -f *.o; rm -f  core  core.*;rm -f findap fdebug 

depclean: clean
	rm -f $(DEPEND)
