![Header](media/header.png)

## 1. INTRODUCTION

In this practice you will learn how to use PyDock [1] to perform docking on a real case from CAPRI experiment [2], target 26 [3]. You will rank 100 docking solutions using pyDock energy, apply experimental data restraints and perform interface prediction based on desolvation energy (Optimal Docking Area [4]) to characterize different properties of the complex. Then, you will have to make your choice and select what you consider to be the best model from our starting pool of docking poses. Finally, RMSD comparison with the real 3-D complex structure will be done to check the results of our *“simulated CAPRI experiment”*.

### 1.1. Download data

Before continuing, clone the data files from this repository in your machine (you will need to install *git* beforehand):

```bash
git clone https://github.com/brianjimenez/pydock_tutorial.git
```

### 1.2. PyDock

We will use here the version 3.0 of pyDock that is already installed in your machine. Test the following in a terminal:

```bash
pydock3
```

If you see an output like this, pyDock is correctly installed:

```bash
PyDock3
  A set of tools for protein-protein docking
  Version 3.5.1

[pyDock3] ERROR: wrong parameters
Usage: pyDock3 dock_name module_name [options]
  pyDock available modules are :
  setup, zdock, ftdock, rotzdock, rotzdock3.0.2, rotftdock, rotftdock2.1, rotref, rotpatchdock, dockser, bindEy, docktet, dockrst
  RMSD, oda, patch, opra, makePDB, makePDBftdock, sipper, show, capriRMSD, randomspin, saxs, resEnergy, dockserContact, dnascore

pyDock3 terminated with error
```

For other installations, you can find further details and instructions in this [web site](https://life.bsc.es/pid/pydock/).


### 1.3. PyDock general syntax

All pyDock jobs are launched as follows:

```bash
pyDock3 dockname module
```

In our example, *dockname* is arbitrarily chosen by the user. In this example, we will use **T26** as *dockname*. The different modules that can be used in pyDock are listed here:

#### Docking:

|    Module        |        Input files       |      Output files                             |
| ---------------- |:------------------------:| :--------------------------------------------:|
|   setup          |      dockname.ini        | dockname\_rec.pdb<br>dockname\_lig.pdb          |
|   ftdock         |      dockname\_rec.pdb<br>dockname\_lig.pdb        | dockname.ftdock       |
|   zdock          |      dockname\_rec.pdb<br>dockname\_lig.pdb        | dockname.zdock        |
|   rotftdock      |      dockname.ftdock       |  dockname.rot  |
|   rotzdock       |      dockname.zdock       |  dockname.rot  |
|   dockser        |      dockname\_rec.pdb<br>dockname\_lig.pdb<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.amber<br>dockname\_rec.pdb.amber<br>      |  dockname.ene  |

#### Additional tools:

|    Module        |        Input files       |      Output files                             |
| ---------------- |:------------------------:| :--------------------------------------------:|
|   dockrst          |      dockname.ini<br>dockname\_rec.pdb<br>dockname\_lig.pdb<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.amber<br>dockname.rot<br>dockname.ene      | dockname.eneRST<br>dockname.rst          |
|   patch         |dockname\_rec.pdb<br>dockname\_lig.pdb<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.H<br>dockname\_rec.pdb.amber<br>dockname.rot<br>dockname.ene      | dockname.recNIP<br>dockname\_rec.pdb.nip<br>dockname\_lig.pdb.nip<br>dockname.ligNIP|
|   oda          | X.pdb<br>Y.pdb      | dockname\_rec.pdb.oda<br>dockname\_rec.oda.ODATAB<br>dockname\_lig.pdb.oda<br>dockname\_lig.oda.ODATAB|


## 2. DOCKING CALCULATIONS

### 2.1. Setup process

The first step before any docking calculation is to generate the pdb files that pyDock will use for docking.

At this point, you must define a receptor and a ligand in your complex. In general, the largest molecule is defined as the receptor and will be kept static, whereas the ligand will be rotated and translated around it.

To begin, you must have in your starting directory:

* [practical/1C5K.pdb](practical/1C5K.pdb)  (TolB or the receptor protein)
* [practical/1OAP.pdb](practical/1OAP.pdb)  (Pal or the ligand protein)

**NOTE: You will find these files at your working directory after `git clone` command.**


In addition, you need an *".ini"* file, which contains the information about the chains to dock from each pdb file, in order to create a new pair of parsed pdb files suitable for pyDock.

Thus, you have to create a text file called **T26.ini** in the `practical` folder and edit it to contain the following (incomplete) information:

```
[receptor]
pdb		= 1C5K.pdb
mol		= A
newmol	= A

[ligand]
pdb		= ?????
mol		= ?????
newmol	= ?????
```

Now we need to complete the *".ini"* file and replace the ***?????*** parts. The `mol` chain name is the original chain name in the input pdb files whereas the `newmol` will be the new chain name in the parsed pdb files suitable for pyDock. 

**Note:** The `newmol` chain names must be different for the receptor and the ligand (so that you can distinguish the chains when docked).

The pdb names in the *“.ini”* file must correspond to the exact names of the pdb files you have in the practical folder (1C5K.pdb, 1c5k.pdb, pdb1c5k.ent.Z, etc...).

#### Remarks:

> * If a pdb does not contain any chain name, use *“-”* in the `mol` field of your *.ini* file.
> * If it contains several copies of the same protein, select only one copy by its chain name.
> * If a protein to dock contains several chains (for example `L` and `H` chains for antibodies) that are relevant for docking, choose different `newmol` names separated by commas.

<br>
Once you have a complete `T26.ini` file, run the pyDock setup writing the following line in your console:

```bash
pydock3 T26 setup
```

This command will create the new PDB files for receptor and ligand:

* T26\_rec.pdb
* T26\_rec.pdb.H (containing hydrogens)
* T26\_rec.pdb.amber (AMBER parameters for each atom in the PDB structure)
* T26\_lig.pdb
* T26\_lig.pdb.H
* T26\_lig.pdb.amber

These files will be suitable for pyDock.


### 2.2. Sampling using Fast Fourier Transform (FFT) methods

PyDock can be applied to score rigid-body docking orientations generated by a variety of methods. We use ZDOCK or FTDock (both are FFT methods) to generate docking positions from T26_rec.pdb and T26_lig.pdb files.








