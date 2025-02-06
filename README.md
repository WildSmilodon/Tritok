# Tritok - a BEM based viscous flow solver

Tritok simulates flow and heat transfer of an incompressible viscous fluid in 3D. Velocity-vorticity formulation of Navier-Stokes eqautions is solved by the Boundary Element Method. It was written by Jure Ravnik, University of Maribor, Faculty of Mechanical Engineering, Smetanova 17, SI-2000, Maribor, Slovenia

Numerical aspects of the algorithm are explained in the following papers:

J. Ravnik, L. Skerget and Z. Zunic: Velocity-vorticity formulation for 3D natural convection in an inclined enclosure by BEM; International Journal of Heat and Mass Transfer, 51, 4517-4527, 2008; http://dx.doi.org/10.1016/j.ijheatmasstransfer.2008.01.018

J. Ravnik, L. Skerget and Z. Zunic: Combined single domain and subdomain BEM for 3D laminar viscous flow Eng. Anal. Bound. Elem., 2008, http://dx.doi.org/10.1016/j.enganabound.2008.06.006

 
### Uses mpich library, compiles under gfortran

```
sudo apt-get update
sudo apt install gfortran make mpich
```


### Installation & test:

```
git clone https://github.com/WildSmilodon/Tritok.git
cd Tritok/src
make
cd ../run
./runAll
```


### Paraview and/or Tecplot are used for results postprocessing

* Install Paraview

```
sudo apt-get update
sudo apt-get install paraview
```

* Open *.vtu files with Paraview to visualise results


### Input files

* ```tri.inp``` - main input file
* ```.bic``` - boundary and initial conditions file
* ```.geo``` - mesh file

### Output files

* ```tri.err``` - error, warnings and messages file (output)
* ```tri.dat``` - results file in ascii Tecplot format (output)
* ```tri.vtu``` - results file in ascii Paraview format (output)
* ```tri.vtk``` - results file in ascii Paraview format (output)
* ```tri.ite``` - solver iteration file (status)
* ```tri.log``` - log file (status)
* ```tri.sta``` - convergence monitoring file (status)
* ```tri.tim``` - time step convergence (status)
* ```tri.tfl``` - temperature flux file (result)
* ```tri.wfl``` - vorticity flux file (result)
* ```tri.int``` - integrals file (input)

### Meshing

GiD (https://www.gidsimulation.com) is used to produce mesh files. Custom addon is needed. The addon and the installation instructions are in the ```docs``` folder.