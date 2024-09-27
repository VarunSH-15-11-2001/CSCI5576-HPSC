# -----------------------------------------------------
# On local Ubuntu or macOS
# -----------------------------------------------------
#
# Be sure to use MPI compiler wrappers.
# Whenever compiling an MPI program, you should use the MPI wrappers:
#
#    C - mpicc
#    C++ - mpiCC, mpicxx, mpic++
#
# These wrappers do all of the dirty work for you of making sure that
# all of the appropriate compiler flags, libraries, include directories,
# library directories, etc. are included when you compile your program.


# -----------------------------------------------------
# Set compiler based on platform
# -----------------------------------------------------

# Default compiler is Intel MPI compiler (Linux)
COMPILER := mpicc

# Check if the system is macOS
ifeq ($(shell uname), Darwin)
    COMPILER := mpic++
endif

# On Ubuntu, change compiler if using a specific user
ifeq ($(strip $(LOGNAME)),scott)
COMPILER := mpic++ 
endif

# -----------------------------------------------------
# Make fd
# -----------------------------------------------------

fd: fd.cpp fd.h plotter.h linear_solver.h
	$(COMPILER) fd.cpp -g -lstdc++ -o fd
