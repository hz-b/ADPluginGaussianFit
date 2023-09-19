# ADPluginGaussianFit

An EPICS areaDetector Plugin to perform 2D Gaussian fit using the optimization problem solver from [Ceres Solver](http://ceres-solver.org/index.html).

Primary Author: Luca Porzio - luca.porzio@helmholtz-berlin.de

## Installation and dependencies

### *Install Ceres Solver*
The ADPluginGaussianFit works with the support of [Ceres Solver](http://ceres-solver.org/index.html) library. It can be build from sources and installed in a selected folder. Please follow the [instructions](http://ceres-solver.org/installation.html) from the official website.


### *Clone Plugin repository*
Clone `ADPluginGaussianFit` repository into your `areaDetector` installation folder, on the same level where `ADCore` is located.

### *Configure areaDetector*
Setup the `areaDetector` installation by adding support to **Ceres Solver**.

In `configure/CONFIG_SITE.local` and add the following lines:
```
WITH_CERES    	= YES
CERES_EXTERNAL	= YES
```
In `configure/CONFIG_SITE.local.<architecture>` and add the following lines:
```
WITH_CERES	= YES
CERES		= /usr/local
#CERES_LIB	= $(CERES)/lib
#CERES_INCLUDE	= $(CERES)/include
```
In `configure/RELEASE_PRODS.local` and add the following lines:
```
# ADPLUGINGAUSSIANFIT is optional plugin
ADPLUGINGAUSSIANFIT=$(AREA_DETECTOR)/ADPluginGaussianFit
```

### *Configure ADCore*
Setup `ADCore` to work with `ADPluginGaussianFit`.

In `commonLibraryMakefile` add the followings lines:
```
ifeq ($(WITH_CERES),YES)
  ifeq ($(CERES_EXTERNAL),NO)
    LIB_LIBS += ceres
  else
    ifdef CERES_LIB
      ceres_DIR        = $(CERES_LIB)
      LIB_LIBS     += ceres
      LIB_LIBS     += glog
      LIB_LIBS     += cholmod
    else
      LIB_SYS_LIBS += ceres
      LIB_SYS_LIBS += glog
      LIB_SYS_LIBS += cholmod
    endif
  endif
endif
```

In `commonDriverMakefile` add the following lines:
```
# ADPLUGINGAUSSIANFIT
ifdef ADPLUGINGAUSSIANFIT
  $(DBD_NAME)_DBD += NDPluginGaussianFit.dbd
  PROD_LIBS	  += NDPluginGaussianFit
  ifdef CERES_LIB
    ceres_DIR +=$(CERES_LIB)
    PROD_LIBS       += ceres glog cholmod
  else
    PROD_SYS_LIBS   += ceres glog cholmod
  endif
endif
```

In `commonPlugins.cmd` add the following lines:
```
NDGaussianFitConfigure("GAUSSIAN_FIT", $(QSIZE), 0, "$(PORT)", 0, 0, 0, 0, 0)
dbLoadRecords("$(ADPLUGINGAUSSIANFIT)/db/NDPluginGaussianFit.template",  "P=$(PREFIX),R=GaussianFit1:, PORT=GAUSSIAN_FIT,ADDR=0,TIMEOUT=1,NDARRAY_PORT=$(PORT)")
set_requestfile_path("$(ADPLUGINGAUSSIANFIT)/db")
```

### *Install Plugin*
Clean and compile in the following order `ADCore`, `ADSupport` and `ADPluginGaussianFit` by running the following command in each directory:
```
make clean uninstall && make -j16
```

### *Do not forget*
In your application, make sure to setup the `ADPLUGINGAUSSIANFIT` variable in the corresponding `envPaths`.