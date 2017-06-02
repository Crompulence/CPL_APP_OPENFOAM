.PHONY: all test clean
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

test:
	@py.test -v ./test
