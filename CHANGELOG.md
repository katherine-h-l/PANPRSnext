# Changes

- Rewrote major portion of the R code
  - Modified files: PANPRS.R PANPRS.cpp
  - Description: Rewrote and renamed many parameters to be more consistent. Modified how some parameters are stored (lists -> matrices) so that they can be passed to C++ as an Armadillo matrix. Parameters are now reorganized and succesfully get passed to the placeholder C++ function.
  - Caveat: The placeholder C++ function does not do anything yet. It just returns one of the input parameters currently as a test for passing by value or passing by reference.

- Introduced Rcpp and RcppArmadillo
  - Modified files: PANPRS.R PANPRS.cpp PANPRS.hpp NAMEPSACE man
  - Description: Armadillo is linear algebra toolkit for C++. Rcpp is interface between R and C++ code. Created a template package and slowly am transferring in old PANPRS code one bit at a time so that it still compiles. Armadillo has support for sparse matrices (instead of using boost). 
  - Caveat: Armadillo does not let us pass by reference (pass pointers to the C++ code), so function definitions and parameters must all be re-written.

- Renamed files and modified variable names accross all files
  - Modified files: all of them
  - Description: snake case for variables, camel_Case for functions. Makes variable names more readable, consistent, and less confusing / distracting. Also renamed files to be more consistent with the variable names.

- Code documentation [WIP]
  - Modified files: gsfPEN.cpp gsfPEN.R
  - Description: Added comments to the code to make it easier to understand what is going on. A lot of the code is redundant and the interface between R and C code can be greatly simplified, but must proceed with caution to ensure that the code still works.

- Changing from C style arrays (using pointers and malloc) to boost arrays [WIP]
  - Modified files: gsfPEN.cpp
  - Description: Changed from C style arrays to boost arrays. We want to leverage the fact that many of the arrays are sparse which requires, so we can save memory. Proof of concept is in place, and we have successfully switched to CPP, and compiled the package using boost. However, the full transition to boost arrays is still in progress.

- Major refactoring of the C code [WIP]
  - Modified files: gsfPEN.cpp
  - Description: Removed (commented out for now) unecessary initialization loops.

- Implemented OpenMP parallelization proof of concept
  - New files: Makevars
  - Description: Proof of concept demonstrating OpenMP can be used within the C code and work with R. Makevars tells R to link to the openMP library when it compiles the C code.

- Switched to C++
  - Modified files: gsfPEN.c gsfPEN.h utility.c utility.h
  - Description: Functions that must be callable from R code are enclosed in extern "C" { ... } blocks.

- Renamed package to PANPRS from SummaryLasso
  - Modified files: NAMESPACE