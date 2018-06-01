DIMAGNO
=======

[![DOI](https://zenodo.org/badge/134246471.svg)](https://zenodo.org/badge/latestdoi/134246471)

![DIMAGNO logo](/docs/logo.png)

DIMAGNO is a two-fluid code for the simulation of DIvergent MAgnetic NOzzles for next generation plasma thrusters (Helicon Plasma Thruster, Electron-Ciclotron Plasma Thruster, Applied-Field MPD thruster, VASIMR, etc).
It uses the method of characteristics to integrate the supersonic
ion equations.

The full version of DIMAGNO was developed since 2009 as part of my MSc and then my PhD theses, and has been used extensively to research the operation and the physics of the plasma flow in a magnetic nozzle. The version contained in this repository is a **refactored** version of that code, created recently with the purpose of making the code open source. As such, the code currently 
only includes the basic functionality, and more features will be included gradually. If you are interested in any particular aspect of the code not yet included in the repository, feel free to contact me.

The full description of the model can be found in *E. Ahedo and M. Merino, "Two-dimensional supersonic plasma acceleration in a magnetic nozzle", Physics of Plasmas 17, 073501 (2010)*.

## Installation

Installation requires simply that you 
[download DIMAGNO](https://github.com/mariomerinomartinez/dimagno/archive/master.zip) 
and add the base directory (the one that contains the `+dimagno` directory) 
to your Matlab path.

### Dependencies

A recent version of Matlab is needed to run the code. 
DIMAGNO has been developed in Matlab R2016b Academic version. 

DIMAGNO depends on other Matlab packages that you can download from this
GitHub account: 
[magnetic_field](https://github.com/mariomerinomartinez/magnetic_field),
[fluid_plasma](https://github.com/mariomerinomartinez/fluid_plasma),
[utilities](https://github.com/mariomerinomartinez/utilities),
[logger](https://github.com/mariomerinomartinez/logger). 
These packages must be installed and in your Matlab path to run DIMAGNO.
 
## Quick usage guide

All input parameters are configured in the `simrc` function included in the package. To use the code, copy this file as a template to the directory where you will run your simulation from, and edit it to your liking. 

After creating a `magnetic_field` and a `fluid_plasma` object using the dependencies listed above, you need to create an initial condition `dimagno.ic`object. You can customize many aspects of these objects and the other parameters in your user `simrc` file, including the name of the folder where the output files will be generated and the type of postprocessing functions that will be run (currently limited to two until I finish refactoring the previous codebase).

Your `simrc` file will be used to overwrite any default values given in the package `simrc` file. You can also override values from the command prompt, given name:value pairs.

To run DIMAGNO, pass the path to your `simrc` file to the main program:

```text
[data,solution] = dimagno.dimagno(path_to_your_simr_file);
```

This will create a number of `.mat` files in the subfolder `fronts` with the raw results of the simulation, and two structures: `data` and `solution`.
The output structure `data` contains all the actual preprocessed inputs used by the program; `solution` contains all the postprocessing results. 

*Unforturnately, the current version of DIMAGNO does not come with detailed documentation.*

### Code structure

The new, refactored DIMAGNO code currently consists of three classes:

* `front`: a class to store the points of the integration front in the method of characteristics.
* `ic`: initial condition class that stores all the interpolation libraries and used in the simulation.
* `straightline`: a helper class for intersecting characteristic lines.

There are also 4 subpackages:

* `exit`: exit conditions to terminate the simulation. The user can select one of them in the input file `simrc`. Currently, only two exit conditions are available; more will be refactored from the old codebase soon.
*  `preprocessor`: contains a function to perform some preprocessing of the input
*  `postprocessor`: functions that the user can select, that will be run after the main simulation loop is completed.
* `solver`: the functions that embody the method of characteristics and enable DIMAGNO to perform a simulation.
 
### Testing

The test suite of the refactored version of DIMAGNO is still quite minimal. 
However, a few tests can be 
found in the `/test` subdirectory. These tests can give you a feeling of how to run DIMAGNO on your own. After adding the package to
your Matlab path, you can run all tests by executing `runtests` from this
subdirectory.

## Contributing

If you have any comments for improvement or 
are interested in contributing to the continued 
development of this or any of my other codes, you can contact us through 
[Mario's website](http://mariomerino.uc3m.es/).

## Acknowledging 

This program is the result of substantial effort. It is released as open
source in the hope that it will be useful to other people. If you find it
useful and/or use it in any of your works, I kindly ask you to acknowledge it
by citing the main article of the Akiles2d model,

> Eduardo Ahedo, Mario Merino, "Two-dimensional supersonic plasma acceleration in a magnetic nozzle", Physics of Plasmas 17, 073501 (2010)  [DOI: 10.1063/1.3442736](https://doi.org/10.1063/1.3442736)

and/or citing the code directly as:

> Mario Merino (2018). DIMAGNO code: Divergent Magnetic Nozzle plasma flow integrator, [DOI: 10.5281/zenodo.1257295](https://doi.org/10.5281/zenodo.1257295)
  
## License

Copyright (c) 2018 Mario Merino. 
The software is released as open source with the [MIT License](LICENSE.md).
