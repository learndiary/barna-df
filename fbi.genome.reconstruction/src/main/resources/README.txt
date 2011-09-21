FLUX CAPACITOR
==============


DESCRIPTION
-----------

The Flux Capacitor is a program to quantify expressed sequence features 
(e.g., genes or transcripts) from transcriptome interrogations by 
short read sequencing.


REQUIREMENTS
------------

* java 1.6 or higher installed

* approx 1Gb of RAM (when SORT_IN_RAM false, see below)

* enough hard disk space for whatever you plan to do


CHANGES
---------------

FluxCapacitor 1.0 RC1

    -fixed issue #11 - all chromosomes from the input annotation are considered now
    -fixed issue #12 and issue #30 - installation problems and issues with native 
    library loading: pre-compiled libraries for Linux-, OSX- (Intel) and Windows-
    platforms are provided for download. In most system configurations libraries 
    are loaded without any special requirement, otherwise they can be set via 
    environment variables.
    -fixed issue#19 - Read identifier formats can now be defined by a flexible 
    expression system.
    -fixed issue #27 - Revised file-IO handling. 
	-fixed issue #63 - Limit RAM-memory usage: From Flux Capacitor v1.0 RC1 on, the 
	program's memory consumption can be made constant by user parameter (SORT_IN_RAM 
	false).

GETTING STARTED
---------------

* unpack

 If you can read this, you obviously managed to unpack the Flux Capacitor :) Welcome

 You will find the executable in the bin/ folder

    flux.bat     Windows

    flux         UNIX clones and Mac OSX


* Run

    You can take a look at the command line parameters and available tools using

     flux --help

    The Flux Capacitor is one of the tools in the Flux Package, further information 
    on each of the tools can be retrieved by the command pattern

     flux --t <toolname> --help

    for instance in the case of the Flux Capacitor

     flux -t capacitor --help

* Parameters

    The Flux Capacitor reads its configuration from a parameter file specified with
    the -p command line parameter (see next section)

    To get a list and descriptions for all available parameters, use

    flux -t capacitor --printParameters

    to create a file with an exhaustive list of parameters and their default values,
    pipe the output to a file

     flux -t capacitor --printParameters > myparameters.par


    NOTE that all file parameters, e.g. the location of the annotation or the mapping file,
    can be specified as either absolute or relative path names; in the latter case the path
    is relative to the location of the parameter file.

* Example (XXX Micha edits here XXX)

    To get a complete sample project including annotations and the genome, download

      http://fluxcapacitor.googlecode.com/files/fluxsimulator_demo_yeast-1.0.tar.gz

	and extract the file

	  tar xzvf fluxsimulator_demo_yeast-1.0.tar.gz

	this will create a demo folder with GTF annotations and the genomic sequences +
	a set of parameters files. Now, to start the simulator do:

	flux -t simulator -p yeast_demo.par

* Memory

    In case you run into out of memory issues--especially when setting parameter 
    SORT_IN_RAM to YES--you can increase the memory size used by the Flux 
    Capacitor using the environment variable FLUX_MEM, for example:

    # this sets the upper limit to 2 gig
    export FLUX_MEM="2G"; flux -t capacitor -p ....


    # this also sets the upper limit to 2 gig
    export FLUX_MEM="2048M"; flux -t capacitor -p ....

* LPSolver libraries

    The Capacitor comes with bundled versions of the LPSolve library for Windows, Linux and Mac OS X.
    In case the bundled libraries are not compatible with your system our you want to use a custom version
    of the library, you have to specify two environment variables

    export LPSOLVER_JNI=/usr/local/lib/liblpsolve55j.so

    export LPSOLVER_LIB=/usr/local/lib/liblpsolve55.so

    to specify teh ABSOLUTE paths to both the shared library (liblpsolve55) and the Java Native Interface
    library (liblpsolve55j). Sources can be downloaded from http://sourceforge.net/projects/lpsolve.
    Note that you have to build both the lpsolve library and the java binding.

* All pages related to the program are reachable from the program homepage

	http://flux.sammeth.net
	
* Especially, have a look at the documentation of the Flux Capacitor at

	http://fluxcapacitor.wikidot.com/capacitor

* If you encounter any bugs, use the bugtracking system at

	http://code.google.com/p/fluxcapacitor/issues/list
	
  First check for any known bugs, if you don't find your issue described, log
  in and create a new issue.

* If you have any questions, discuss them via the forum 
	
	http://fluxcapacitor.wikidot.com/forum:start
	
or leave me an email.

The FluxDev-Team

Thasso and Micha
