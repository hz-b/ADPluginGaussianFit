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

ifdef CERES_INCLUDE
    USR_INCLUDES += -I$(CERES_INCLUDE)
endif
ifdef CERES_LIB
    NDPluginGaussianFit_DIR += $(CERES_LIB)
endif

USR_INCLUDES += -I/usr/include/eigen3
USR_INCLUDES += -I/usr/include/glog

include $(AREA_DETECTOR)/ADCore/ADApp/commonLibraryMakefile

include $(TOP)/configure/RULES