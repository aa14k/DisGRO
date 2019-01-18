Distance-guided Sequential chain-Growth Monte Carlo(DiSGro)

DiSGro is a fast protein loop sampling and structure prediction tool
based on distance-guided sequential chain-growth Monte Carlo method.
It can efficiently generate high quality near-native loop
conformations which can be used as candidates for further refinement.

If you use DiSGro, cite the following reference:

Ke Tang, Jinfeng Zhang, and Jie Liang, (2014) Fast protein loop
sampling and structure prediction using distance-guided sequential
chain-growth Monte Carlo method. PLoS Computational Biology. In press.


1. Installation DiSGro:

    Untar DiSGro.tar.gz 

      tar -xvzf DiSGro.tar.gz

    Compile with
      make

2. Data Preparation

   Input files must be in the PDB format. The atoms of absent loop
   regions need to be labeled as "H" atoms, whose coordinates are assigned to
   0.000 for each axis.  

   See the provided example - 1ctqa.pdb.

3. Running DiSGro: 

./disgro -mode [Mode] -f [Input PDB] -n [Input # of trials] -nds [# of sampling states] -start [Start position] -end [End position] -confkeep [# of retained conformations] -ellip [REsidue-residue Distance Cutoff and ELLipsoid criterion (Redcell)] -nscc [# of side chain states] -pdbout [# of output models]


Description:

  -mode       mode of the computation, currently mode: smc.
              smc: loop modeling using sequential chain-growth Monte Carlo
  -f          coordinate file of the protein (.pdb)
  -n          number of loop conformations trying to generate
  -nds        number of states sampling from an empirically derived backbone dihedral angle distribution. It corresponds to n in PLOS Computational Biology paper, and it is a subset  of m, which is the number of states sampled from our distance distributions
  -start      the start position for growth mode.
  -end        the end position for growth mode.
  -confkeep   number of retained conformations for adding side-chains
  -eval       energy evalutaion of conformations
  -ellip      turn on or turn off REsidue-residue Distance Cutoff and ELLipsoid criterion (Redcell)
  -nscc       number of side chain states
  -pdbout     number of output conformations
  -kcluster   k-mean clustering


Example:

./disgro -mode smc -f 1ctqa.pdb -n 5000 -nds 32 -start 26 -end 37 -eval -confkeep 1000 -ellip -nscc 20 -pdbout 100

4. Getting RMSD:

After running Disgro, the RMSD for the conformations will be stored in a file called RMSD.pdb. The first column of the file is the conformation number while the second column is the RMSD(original protien, Disgro conformation). 


Contact
ktang6@uic.edu
aa14k@my.fsu.edu
