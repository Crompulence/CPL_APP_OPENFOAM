.PHONY: all test clean
all:
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam
	
clean:
	@wclean src/CPLSocketFOAM
	@wclean src/solvers/CPLIcoFoam
	rm -rf bin
	rm -rf lib

test:
	@py.test -v ./test

patch-scotch:                                                                                                                                                                                                                                    
	patch $(FOAM_SRC)/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C ./config/ptscotchDecomp.patch
