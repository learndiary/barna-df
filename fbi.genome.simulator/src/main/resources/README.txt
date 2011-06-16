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


GETTING STARTED
---------------

* unpack and run the batch file in the /bin directory that corresponds to 

your platform:

flux.bat     Windows

flux         UNIX clones and Mac OS X


* For having some sample data, download the corresponding demo data

    http://code.google.com/p/fluxcapacitor/issues/list

    For example, get :
    http://fluxcapacitor.googlecode.com/files/fluxsimulator_demo_yeast-1.0.tar.gz

	and extract the file

	tar xzvf fluxsimulator_demo_yeast-1.0.tar.gz

	this will create a demo folder with GTF annotations and teh genomic sequences +
	a set of parameters files. Now, to start the simulator (assuming that flux is in
	your path) do:

	flux -t simulator -p yeast_demo.par


* All pages related to the program are reachable from the program homepage

	http://flux.sammeth.net
	
* Especially, have a look at the documentation of the FLUX SIMULATOR at

	http://fluxcapacitor.wikidot.com/simulator
	

* If you encounter any bugs, use the bugtracking system at

	http://code.google.com/p/fluxcapacitor/downloads/list
	
First check for any known bugs, if you don't find your issue described, log 

in and create a new issue.


* If you have any questions, discuss them via the forum 
	
	http://fluxcapacitor.wikidot.com/forum:start
	
or leave me an email.


The FluxDev-Team