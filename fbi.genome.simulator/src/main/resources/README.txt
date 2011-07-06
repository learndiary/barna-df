FLUX SIMULATOR
==============


DESCRIPTION
-----------

The Flux Simulator is a program to simulate in silico RNAseq experiments as
carried out adopting new sequencing technologies.


REQUIREMENTS
------------

* java 1.6 or higher installed

* approx 1Gb of RAM

* enough hard disk space for whatever you plan to do


CHANGES
---------------

FluxSimulator RC2

    -fixed issue #55 and allow to disable TSS_MEAN by setting it to NaN
    -fixed issue #57 - the GTF sorter
    -added an amplification step and removed RT_GC_LO parameter. Now, GC content
     is weighted in amplification. The amplification consists of three parameters:
      * PCR_ROUNDS - the number of rounds (set this to 0 to disable amplification)
      * GC_MEAN and GC_SD - mean and standard deviation for the GC content distribution


GETTING STARTED
---------------

* unpack

 If you can read this, you obviously managed to unpack the FluxSimulator :) Welcome
 (not necessarily vice versa, sufficient condition)

 You will find the executable in the bin/ folder

    flux.bat     Windows

    flux         UNIX clones and Mac OS X


* Run

    You can take a look at the command line parameters and available tools using

     flux --help

    To start a specific tool and see its parameters, use

     flux --t <toolname> --help

    for example

     flux -t simulator --help

* Parameters

    The FluxSimulator reads its configuration from a parameter file specified with
    the -p command line parameter (see next section)

    To get a list and descriptions for all available parameters, use

    flux -t simulator --printParameters

    to create a file with an exhaustive list of parameters and their default values,
    pipe the output to a file

     flux -t simulator --printParameters > myparameters.par


    NOTE that all file parameters, e.g. the location of the genome or the .gtf annotation file,
    can be specified as either absolute or relative path names; in the latter case the path
    is relative to the location of the parameter file.
    Also, in contrast to previous versions, all result files (i.e. .pro, .lib, .bed, .fasta)
    are created in the same directory and with the same name as the parameter file, for instance:

    foo.par

    creates the default results as

    foo.pro, foo.lib, and foo.bed

    You can still explicitly set the output file names in the parameter file with corresponding values
    for e.g. PRO_FILE_NAME.

* Example

    To get a complete sample project including annotations and the genome, download

      http://fluxcapacitor.googlecode.com/files/fluxsimulator_demo_yeast-1.0.tar.gz

	and extract the file

	  tar xzvf fluxsimulator_demo_yeast-1.0.tar.gz

	this will create a demo folder with GTF annotations and the genomic sequences +
	a set of parameters files. Now, to start the simulator do:

	flux -t simulator -p yeast_demo.par

* Memory

    In case you run into out of memory issues, you can increase the memory size used by
    the simulator using the environment variable FLUX_MEM, for example:

    # this sets the upper limit to 2 gig
    export FLUX_MEM="2G"; flux -t simulator -p ....


    # this also sets the upper limit to 2 gig
    export FLUX_MEM="2048M"; flux -t simulator -p ....

* All pages related to the program are reachable from the program homepage

	http://flux.sammeth.net
	
* Especially, have a look at the documentation of the FLUX SIMULATOR at

	http://fluxcapacitor.wikidot.com/simulator

* If you encounter any bugs, use the bugtracking system at

	http://code.google.com/p/fluxcapacitor/issues/list
	
  First check for any known bugs, if you don't find your issue described, log
  in and create a new issue.

* If you have any questions, discuss them via the forum 
	
	http://fluxcapacitor.wikidot.com/forum:start
	
or leave me an email.

The FluxDev-Team

Micha and Thasso
