build: 
	g++ ../src/weightedMaximumStableSolver.cpp ../src/graphNO.cpp ../src/solver_to_csv.cpp ../test/testWeightedMaximumStableSolver.cpp

test: 
	./a.out 

create: 
	python3 ../script/generation.py

buildMac:
	clang++ -Wall -std=c++11 ../src/weightedMaximumStableSolver.cpp ../src/graphNO.cpp ../src/solver_to_csv.cpp ../test/testWeightedMaximumStableSolver.cpp

