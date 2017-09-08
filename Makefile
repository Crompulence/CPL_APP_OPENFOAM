OpenFOAM_DIR=`cat CODE_INST_DIR`
OpenFOAM_SRC_DIR=$(OpenFOAM_SRC_DIR)/src
OpenFOAM_ETC_DIR=$(OpenFOAM_ETC_DIR)/etc

.PHONY: all test clean clean-test
all:
	@wmake libso src/CPLPstream
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam
	
patch-openfoam:
	cp ./config/pref.sh $(OpenFOAM_ETC_DIR)/config/

clean:
	@wclean src/CPLSocketFOAM
	@wclean src/solvers/CPLIcoFoam
	@wclean src/CPLPstream
	rm -rf bin
	rm -rf lib

clean-test:
	cd test/stressC-P/debug && ./clean.sh
	cd test/velocityP-C/debug && ./clean.sh

test:
	@py.test2 -v ./test
