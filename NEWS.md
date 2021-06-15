# FMM 0.3.0
Enhancements:
* Optimization of the estimation procedure for the monocomponent, multicomponent and restricted FMM models. An embedded parallelized procedure is available for all the models for faster estimations.
* Code has been refined: it is now more object-oriented, code duplicity has been reduced, unnecessary orders have been suppressed, API usability has been improved, and the documentation is now automatically created with roxygen2.

New features:
* In addition to the estimation by blocks previously available for the restricted model, an exact implementation has been added for more accurate fits in expense of more computational time.

# FMM 0.2.0

New features:
* Optimized monocomponent and multicomponent model computation.

# FMM 0.1.2

Enhancements:
* Added a slot in the FMM object with information on the number of iteraction of the fitting process.
* Renamed the FMMPeaks function to getFMMPeaks.
* Renamed the generate_FMM function to generateFMM. Also parameters a, b and w to alpha, beta and omega, respectively.

New features:
* New show() method to print basic information on an FMM object.

# FMM 0.1.1

Fixed the issues pointed out in the review of the CRAN submission:

* Modified the beginning of the description.
* Added details about internal functions in the documentation.
* Replaced the use of cat() by warning().
* options(warn=-1) is no longer used.
* No example in the documentation of the functions uses more than two processor cores.

# FMM 0.1.0

* Submitted version 0.1.0 to CRAN.
