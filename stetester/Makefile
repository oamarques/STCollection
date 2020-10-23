###########################################################################
#                                                                         #
# Usage:                                                                  #
#                                                                         #
#    make [target]                                                        #
#                                                                         #
# where target can be                                                     #
#                                                                         #
#    dstetester.x : double precision version of stetester                 #
#    sstetester.x : single precision version of stetester                 #
#    stetesterlib : creates the stetester library                         #
#                                                                         #
# By default all the three targets are created.                           #
#                                                                         #
###########################################################################

################
include make.inc
################

.SUFFIXES:
.SUFFIXES: .f90 .f .F .o

all: stetesterlib dstetester.x sstetester.x 

# Object files for stetester ##############################################

dstetester_lib = \
	SRC/dstedefinitions.o \
	SRC/mapddata.o \
	SRC/dsteinterfaces1.o \
	SRC/dsteinterfaces2.o \
	SRC/dsteinterfaces3.o \
	SRC/dbdchkrslt.o \
	SRC/dbdprnrslt.o \
	SRC/dstebndgap.o \
	SRC/dstechkrslt.o \
	SRC/dstedefbmtrx.o \
	SRC/dstedeftmtrx.o \
	SRC/dstedefvals.o \
	SRC/dstegetdata.o \
	SRC/dstematgen.o \
	SRC/dstemccallst.o \
	SRC/dstemcdump.o \
	SRC/dstemcglued.o \
	SRC/dstemcmatrix.o \
	SRC/dstemcmatrixf.o \
	SRC/dstemcvalsi.o \
	SRC/dstemcvalsv.o \
	SRC/dstemcvalues.o \
	SRC/dstemcvaluesf.o \
	SRC/dsteprnrslt.o \
	SRC/dsterndiliu.o \
	SRC/dsterndvlvu.o \
	SRC/getdreal.o
   
sstetester_lib = \
	SRC/sstedefinitions.o \
	SRC/mapsdata.o \
	SRC/ssteinterfaces1.o \
	SRC/ssteinterfaces2.o \
	SRC/ssteinterfaces3.o \
	SRC/sbdchkrslt.o \
	SRC/sbdprnrslt.o \
	SRC/sstebndgap.o \
	SRC/sstechkrslt.o \
	SRC/sstedefbmtrx.o \
	SRC/sstedeftmtrx.o \
	SRC/sstedefvals.o \
	SRC/sstegetdata.o \
	SRC/sstematgen.o \
	SRC/sstemccallst.o \
	SRC/sstemcdump.o \
	SRC/sstemcglued.o \
	SRC/sstemcmatrix.o \
	SRC/sstemcmatrixf.o \
	SRC/sstemcvalsi.o \
	SRC/sstemcvalsv.o \
	SRC/sstemcvalues.o \
	SRC/sstemcvaluesf.o \
	SRC/ssteprnrslt.o \
	SRC/ssterndiliu.o \
	SRC/ssterndvlvu.o \
	SRC/getsreal.o

gstetester_lib = \
	SRC/gstedefinitions.o \
	SRC/gsteinterfaces1.o \
	SRC/gsteinterfaces2.o \
	SRC/getintgr.o \
	SRC/getrecord.o \
	SRC/getstrng.o \
	SRC/handler.o \
	SRC/listlength.o \
	SRC/parser.o \
	SRC/parserlist.o

dstetester_obj = SRC/dstetester.o SRC/dsteflops.o SRC/dstetimer.o
sstetester_obj = SRC/sstetester.o SRC/ssteflops.o SRC/sstetimer.o

# Target _stetester.x requires LAPACK #####################################

dstetester.x: $(dstetester_obj) SRC/dsteruntests.o
	$(FC) $(FCLDFLGS) -o $@ $^ $(STETESTERLIB) \
	$(LAPACKLIB) $(TMGLIB) $(BLASLIB)

sstetester.x: $(sstetester_obj) SRC/ssteruntests.o
	$(FC) $(FCLDFLGS) -o $@ $^ $(STETESTERLIB) \
	$(LAPACKLIB) $(TMGLIB) $(BLASLIB)

# Target stetesterlib #####################################################

stetesterlib: $(gstetester_lib) $(dstetester_lib) $(sstetester_lib)
	$(ARCH) $(ARCHFLAGS) $(STETESTERLIB) $^
	$(RANLIB) $(STETESTERLIB)

###########################################################################

clean:
	-@rm -f SRC/*.o *.log *.m *_dump.*

cleanall:
	-@rm -f *.x LIB/*.a LIB/*.mod; make clean

SRC/dstetimer.o: SRC/dstetimer_$(TIMER).F
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) -c -o $@ $<
SRC/sstetimer.o: SRC/sstetimer_$(TIMER).F
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) -c -o $@ $<
SRC/dsteflops.o: SRC/dsteflops_$(FLOPS).F
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) -c -o $@ $<
SRC/ssteflops.o: SRC/ssteflops_$(FLOPS).F
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) -c -o $@ $<

.f90.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) -c -o $@ $< 
.F.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) -c -o $@ $<
