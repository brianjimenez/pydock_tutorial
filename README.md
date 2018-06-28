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
$ pydock3
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
$ pyDock3 dockname module
```

In our example, *dockname* is arbitrarily chosen by the user. In this example, we will use **T26** as *dockname*. The different modules that can be used in pyDock are listed here:

