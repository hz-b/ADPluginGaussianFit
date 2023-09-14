TOP=..
include $(TOP)/configure/CONFIG

DIRS := $(DIRS) $(filter-out $(DIRS), $(wildcard *db*))
# DIRS := $(DIRS) $(filter-out $(DIRS), $(wildcard *Db*))
# DIRS := $(DIRS) $(filter-out $(DIRS), $(wildcard *op*))

DIRS += src

include $(TOP)/configure/RULES_DIRS