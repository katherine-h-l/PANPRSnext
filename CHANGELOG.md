# Changes

- Renamed package to PANPRS from SummaryLasso
  - Modified files: NAMESPACE
  
- Switched to C++
  - Modified files: gsfPEN.c gsfPEN.h utility.c utility.h
  - Description: Functions that must be callable from R code are enclosed in extern "C" { ... } blocks.
  
- Implemented OpenMP parallelization proof of concept
  - New files: Makevars
  - Description: Proof of concept demonstrating OpenMP can be used within the C code and work with R. Makevars tells R to link to the openMP library when it compiles the C code.

- Major refactoring of the C code [WIP]
  - Modified files: gsfPEN.cpp
  - Description: Removed (commented out for now) unecessary initialization loops.

- Changing from C style arrays (using pointers and malloc) to boost arrays [WIP]
  - Modified files: gsfPEN.cpp
  - Description: Changed from C style arrays to boost arrays. We want to leverage the fact that many of the arrays are sparse which requires, so we can save memory. Proof of concept is in place, and we have successfully switched to CPP, and compiled the package using boost. However, the full transition to boost arrays is still in progress.

- Code documentation [WIP]
  - Modified files: gsfPEN.cpp gsfPEN.R
  - Description: Added comments to the code to make it easier to understand what is going on. A lot of the code is redundant and the interface between R and C code can be greatly simplified, but must proceed with caution to ensure that the code still works.