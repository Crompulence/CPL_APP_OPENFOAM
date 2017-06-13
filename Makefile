.PHONY: all test clean clean-test
all:
	@wmake libso src/CPLPstream
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam

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
	@py.test -v ./test
