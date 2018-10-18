## Installation

Prereqesuites:
The package only needs pacakges available in the `cvmfs` and the on GitHub publicly available DooSoftware pacakges:

 * DooCore: https://github.com/e5-tu-do/DooCore
 * DooFit: https://github.com/e5-tu-do/DooFit
 * DooSelection: https://github.com/e5-tu-do/DooSelection

To install the package then you need first to check out the repository with all its submodules (in the example the access via ssh is shown):

```
git clone ssh://git@git.e5.physik.tu-dortmund.de:10022/alex.birnkraut/CPTimeFit_BD2Dpi.git
cd CPTimeFit_BD2Dpi
git submodule update --init --recursive
```

For the further installation it is recommended to create a build directory to store all cmake/make files which are created. Then the commands to build the package can be executed:

```
mkdir build
cd build
cmake ..
make main
```

This creates an exeuctable at `bin/main` which can be configured via a configfile (see examples in the `config_Bd2DpiMC` and `config_Bd2DpiData` directories).

## Usage

To execute the Fitter one can e.g. execute the following command (assuming that you're still in the `build` directory):

```
./bin/main -C ConfigFile
```
