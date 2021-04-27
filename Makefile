# This file is run by the command `make`
# with `make` the file starts with the first target

# Directory variables
DIR		= ./
IDIR	= $(DIR)include/
MDIR	= $(DIR)modules/
TDIR	= $(DIR)tests/

# variables of flags
CC = gcc
OPTIMIZATION = -O2
# -Wall gcc warnings during compilation
CFLAGS = -Wall -lm $(OPTIMIZATION)

# variables with executable files
TEST_OBS = test_observables
TEST_FULLMAG = testfullmag
TEST_GEO = testgeometry
TEST_CON = testconsistency
MAIN = main

# variable for the modules
# $(wildcar *) is the safe version of * and means all files in the modules
MODULES = $(wildcard $(MDIR)*.c)

# variables for including
INCLUDE = -I $(IDIR)

# this is the first target all with the depency OBJECT_FILES
# it will look for the depency before running the below command(s)
all:
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(MAIN).c -o $(MAIN).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_LIN).c -o $(TEST_OBS).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_ANAL).c -o $(TEST_FULLMAG).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_UNI).c -o $(TEST_GEO).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_GAUSS).c -o $(TEST_CON).exe

# maybe its better to compile and link seperately
# all: $(OBJECT_FILES)
# 		echo: "Linking: $@ ($(CC))"
# 		$(CC) -o $(TEST_LIN)

# $(OBJECT_FILES):
