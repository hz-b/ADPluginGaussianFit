# ADPluginGaussianFit

An EPICS areaDetector Plugin to perform 2D Gaussian fit using the optimization problem solver from [Ceres Solver](http://ceres-solver.org/index.html).

Primary Author: Luca Porzio - luca.porzio@helmholtz-berlin.de

## Installation and dependencies

The ADPluginGaussianFit works with the support of [Ceres Solver](http://ceres-solver.org/index.html) library. It can be build from sources and installed in a selected folder. Please follow the [instructions](http://ceres-solver.org/installation.html) from the official website.


After the installation, enter the `ADCore` directory, and open `ADApp/commonDriverMakefile`. Ensure that the following is added and uncommented after `ADPluginEdge` or `ADPluginBar`:

```makefile
ifdef ADPLUGINGAUSSIANFIT
  $(DBD_NAME)_DBD += NDPluginGaussianFit.dbd
  PROD_LIBS	  += NDPluginGaussianFit
  ifdef CERES_SOLVER_LIB
    ceres_DIR +=$(CERES_SOLVER_LIB)
    PROD_LIBS       += ceres
  else
    PROD_SYS_LIBS   += ceres
  endif
endif
```

