OpenFOAM_DIR=`cat CODE_INST_DIR`
OpenFOAM_SRC_DIR=$(OpenFOAM_SRC_DIR)/src
OpenFOAM_ETC_DIR=$(OpenFOAM_ETC_DIR)/etc

.PHONY: all test clean clean-test
all:
	@wmake libso src/CPLPstream
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam
	@wmake src/solvers/CPLporousIcoFoam

sedifoam:
	@wmake libso src/CPLPstream
	@wmake libso src/CPLSocketFOAM
	@wmakeLnInclude src/solvers/CPLSediFoam/dragModels
	@wmake libso src/solvers/CPLSediFoam/dragModels
	@wmakeLnInclude src/solvers/CPLSediFoam/chPressureGrad
	@wmake libso src/solvers/CPLSediFoam/chPressureGrad
	@wmakeLnInclude src/solvers/CPLSediFoam/lammpsFoamTurbulenceModels
	@wmake libso src/solvers/CPLSediFoam/lammpsFoamTurbulenceModels
	@wmake src/solvers/CPLSediFoam
	
patch-openfoam:
	cp ./config/pref.sh $(OpenFOAM_ETC_DIR)/config/

clean:
	@wclean src/CPLSocketFOAM
	@wclean src/solvers/CPLIcoFoam
	@wclean src/solvers/CPLporousIcoFoam
	@wclean src/CPLPstream
	@wclean src/CPLSediFoam
	rm -rf bin
	rm -rf lib

clean-test:
	cd test/stressC-P/debug && ./clean.sh
	cd test/velocityP-C/debug && ./clean.sh

test:
	@py.test2 -v ./test

#.PHONY: all test clean clean-test
#all: background CPLIcoFOAM CPLporousIcoFoam
#	@echo "Building everything"

#background: 
#	@wmake libso src/CPLPstream
#	@wmake libso src/CPLSocketFOAM

#CPLIcoFOAM:
#	@wmake src/solvers/CPLIcoFoam

#CPLporousIcoFoam:
#	@wmake src/solvers/CPLporousIcoFoam
#	
#patch-openfoam:
#	cp ./config/pref.sh $(OpenFOAM_ETC_DIR)/config/

#clean:
#	@wclean src/CPLSocketFOAM
#	@wclean src/solvers/CPLIcoFoam
#	@wclean src/CPLPstream
#	rm -rf bin
#	rm -rf lib

#clean-test:
#	cd test/stressC-P/debug && ./clean.sh
#	cd test/velocityP-C/debug && ./clean.sh

#test:
#	@py.test -v ./test

#patch-scotch:                                                                                                                                                                                                                                    
#	patch $(FOAM_SRC)/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C ./config/ptscotchDecomp.patch

