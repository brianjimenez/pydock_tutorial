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
pdb     = 1C5K.pdb
mol     = A
newmol  = A

[ligand]
pdb     = ?????
mol     = ?????
newmol  = ?????
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

PyDock can be applied to score rigid-body docking orientations generated by a variety of methods. 

We use ZDOCK or FTDock (both are FFT-based methods) to generate docking positions from T26_rec.pdb and T26_lig.pdb files.

Usually, we will use the `ftdock` or `zdock` modules to run the sampling using FTDock or ZDOCK software:

```bash
pydock3 T26 ftdock
```

or

```bash
pydock3 T26 zdock
```

But in this practical we will skip this step as it may require a long time to compute (from several minutes to hours depending on the CPU cores available).

So now choose between using the results from ZDOCK or FTDock and copy the need files to your `practical` folder.

#### Using FTDock results

```
cd practical
cp ../T26/dockings/ftdock/T26.ftdock .
```

#### Using ZDOCK results

```
cd practical
cp ../T26/dockings/zdock/T26.zdock .
```

#### 2.2.1. From FFT to Rot

Then, for each conformation, we need to transform the output data from ZDOCK (T26.zdock) and FTDock (T26.ftdock) (in which each solution is represented by the cartesian position of ligand and the rotation based on Euler angles) to the rotation and translation matrix that transforms the original ligand into the conformations generated by FTDock or ZDOCK. 

This is done by using the following command for ZDOCK:

```bash
pydock3 T26 rotzdock
```

or this one for FTDock:

```bash
pydock3 T26 rotftdock
```

This calculation is quite fast and it will create a T26.rot file containing the transformation matrices above mentioned for all docking poses. 

**IMPORTANT: Because of time limitations, in this practical we will only proceed with a subset of docking solutions. For that, you need to edit the T26.rot file to keep just 100 conformations (you may choose randomly 100 lines from the file).**

Before editing T26.rot file, make a copy of the original T26.rot file to save the entire set of docking solutions:

```bash
cp T26.rot T26.rot.orig
```

**Tip:** you can select *N* lines from a file using `sed`:

```bash
sed -n -e '10,110p' T26.rot.orig > T26.rot
```

### 2.3. Scoring using the pyDock energy function

Next stage is to use pyDock energy function to score and rank all positions by running dockser module with the following command:

```bash
pydock3 T26 dockser  > dockser.log  &
```

**WARNING: Be sure you created a T26.rot file with just 100 conformations, otherwise this step would last several hours.**

When dockser finishes, take a look to the output file called `T26.ene` that will look like the following example, with different values:

```
      Conf(1)     Ele(2)     Desolv (3)   VDW(4)     Total(5)        RANK(6)
----------------------------------------------------------------------------------
        8726     -28.979      -9.712     130.111      -38.691             1
        4538     -28.001      -8.980      38.482      -36.981             2
        6446     -29.716      -4.215      96.438      -33.931             3
        1590     -32.394       0.109      28.699      -32.285             4

```

* (1)	**Conf:** conformation number of the docking pose (same as that in the rot file, last column)
* (2)	**Ele:** Electrostatic energy component
* (3)	**Delsov:** Desolvation energy component
* (4)	**VDW:** Van der Waals energy component (term weighted to 0.1 by default)
* (5)	**Total:** Total binding energy (Ele + Desolv + 0.1*VDW)
* (6)	**RANK:** conformation rank according to its computed binding energy







