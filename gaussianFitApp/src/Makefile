TOP=../..
include $(TOP)/configure/CONFIG

# The following gets rid of the -fno-implicit-templates flag on vxWorks, 
# so we get automatic template instantiation.
# This is what we want for the source files that define template functions.
ifeq (vxWorks,$(findstring vxWorks, $(T_A)))
    CODE_CXXFLAGS=
endif

DBD += NDPluginGaussianFit.dbd

INC += NDPluginGaussianFit.h

LIBRARY_IOC += NDPluginGaussianFit

NDPluginGaussianFit_SRCS += NDPluginGaussianFit.cpp

ifdef CERES_SOLVER_INCLUDE
    USR_INCLUDES += -I$(CERES_SOLVER_INCLUDE)
endif
ifdef CERES_SOLVER_LIB
    NDPluginGaussianFit_DIR += $(CERES_SOLVER_LIB)
endif

include $(AREA_DETECTOR)/ADCore/ADApp/commonLibraryMakefile

include $(TOP)/configure/RULES