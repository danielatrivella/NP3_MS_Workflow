TRGTS = Common_lib MLlib_lib PepNoVo_lib MsCluster_lib

all: $(TRGTS)

Common_lib:
	cd ./src/Common ; make

MLlib_lib:
	cd ./src/MLlib ; make

PepNoVo_lib:
	cd ./src/PepNovo ; make

MsCluster_lib:
	cd ./src/MsCluster ; make ; cp NP3_MSCluster_bin ../..

Models:
	cd ./Models ; ./dos2unix *.* ; dos2unix */*.*

clean:
	cd ./src/Common ; make clean
	cd ./src/MLlib ; make clean
	cd ./src/PepNovo ; make clean
	cd ./src/MsCluster ; make clean

