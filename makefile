disgro: main.o atom.o residue.o structure.o potential.o smc.o rotamer.o reprst.o util.o cal_energy.o matrix.o kmean_clustering.o sample_states.o 
	g++ -o disgro main.o atom.o residue.o smc.o structure.o potential.o rotamer.o reprst.o util.o cal_energy.o matrix.o sample_states.o kmean_clustering.o -lm -w
disgrodb: main.o atom.o residue.o structure.o potential.o rotamer.o reprst.o util.o cal_energy.o mc_move.o
	g++ -g -o disgrodb main.o atom.o residue.o nns.o structure.o potential.o rotamer.o reprst.o util.o cal_energy.o mc_move.o -lm
main.o: main.cpp atom.h residue.h structure.h rotamer.h potential.h reprst.h util.h cal_energy.h smc.h
	g++ -c -O3 main.cpp 

smc.o: smc.cpp smc.h structure.h residue.h atom.h potential.h util.h
	g++ -c -O3 smc.cpp
cal_energy.o: cal_energy.h cal_energy.cpp structure.h util.h potential.h rotamer.h residue.h atom.h
	g++ -c -O3 cal_energy.cpp
sample_states.o: sample_states.cpp sample_states.h atom.h residue.h structure.h util.h rotamer.h reprst.h
	g++ -c -O3 sample_states.cpp
atom.o: atom.h atom.cpp
	g++ -c -O3 atom.cpp
residue.o: residue.cpp residue.h atom.h util.h
	g++ -c -O3 residue.cpp
structure.o: structure.cpp structure.h atom.h residue.h util.h reprst.h sample_states.h cal_energy.h
	g++ -c -O3 structure.cpp
potential.o: potential.h potential.cpp atom.h residue.h structure.h rotamer.h util.h cal_energy.h
	g++ -c -O3 potential.cpp
rotamer.o: rotamer.cpp rotamer.h atom.h residue.h
	g++ -c -O3 rotamer.cpp
reprst.o: reprst.cpp reprst.h util.h atom.h residue.h
	g++ -c -O3 reprst.cpp
util.o: util.h util.cpp atom.h structure.h residue.h matrix.h
	g++ -c -O3 util.cpp
matrix.o: matrix.h matrix.cpp atom.h
	g++ -c -O3 matrix.cpp
kmean_clustering.o: kmean_clustering.cpp
	g++ -c -O3 kmean_clustering.cpp
depend:
	gcc -E -MM *.cpp > .depend
clean: 
	rm -f *.o *~
	
