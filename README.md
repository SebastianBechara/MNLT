# MNLT 
This is a repository for the Mariusz Nowacki Lab Tool-box. 

It comprises some common workflows bundled together in a python program. 

## Install
This tool box require a bunch of other programs. An installation guide is provided at the begining of the documentation. Once all dependencies are installed you can invoke the `check` scripts to see if everything works properly. `check_I` will just try to find the dependencies and `check_II` will run the full program using a test set. This will take quite a while, also ich will use 6 threads. If you are working on an old machine you may want to check the required programs by hand.

## Usage
A full walkthrough of the different modules is provided in the documentation. 

In short: You can run MNLT using this command
```
./MNLT 
```
Output:
```
############################################################################
A wrapper script for multiple sub-scripts. Contains the following scripts:
  *  IRScalc.py :: Calculates IRS
  *  EWMA_FreqPLOTS.py :: Generates EWMA and nucleotide frequency plots
  *  afterParTIES.py :: Generates Estiennes correlation matrices
  *  Clustering.py :: Alternative for Estiennes correlation matrices using
                     dimensional reduction
  *  sRNAHist.py :: Generates the typical sRNA barplots
  *  ParTIES_MILORD :: MILORD submodule of the software package for
                       paramecium tools
  *  afterMILORD.py :: Generates alternative / cryptic excision plots
############################################################################
Usage:
    MNLT.py <Module> <Module arguments>

Modules:
  IRS_calc      Run IRS calculation
  EWMA_Freq     Run EWMA and Frequency calculation
  afterParties  Run Estiennes correlation analysis
  Clustering    Run clustering using dimension reduction
  sRNA_hist     Create sRNA histogrammes
  MILORD        Run ParTIES submodule MILORD for excision events
  afterMILORD   Create cryotic/alternate excision plots

Optional:
  -h, --help    Show this help message and exit
```
Each submodule is callable like so: `./MNLT [SUBMODULE]`
