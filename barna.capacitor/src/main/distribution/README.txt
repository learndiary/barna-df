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

FluxCapacitor 1.0.1
    - BARNA-126 - added minimal set of command line options

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

    flux-capacitor.bat     Windows

    flux-capacitor         UNIX clones and Mac OSX


* Run

    You can take a look at the command line parameters and available tools using

     flux-capacitor --help

    The Flux Capacitor is one of the tools in the Flux Package, further information 
    on each of the tools can be retrieved by the command pattern

     flux-capacitor --t <toolname> --help


* Parameters

    The Flux Capacitor reads its configuration from a parameter file specified with
    the -p command line parameter (see next section)

    To get a list and descriptions for all available parameters, use

    flux-capacitor --printParameters

    to create a file with an exhaustive list of parameters and their default values,
    pipe the output to a file

     flux-capacitor --printParameters > myparameters.par


    NOTE that all file parameters, e.g. the location of the annotation or the mapping file,
    can be specified as either absolute or relative path names; in the latter case the path
    is relative to the location of the parameter file.

* Example

    To get a complete sample project including annotations and the genome, download

      http://fluxcapacitor.googlecode.com/files/fluxcapacitor_demo_mouse-1.0.tar.gz

	and extract the file

	  tar xzvf fluxcapacitor_demo_mouse-1.0.tar.gz

	this will create a demo folder with GTF annotation, the mapped reads' file and
	a parameter file. Now, to start the capacitor do:

	flux-capacitor -p mouse_demo.par

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

LICENSE
---------------
The Flux Capacitor and all Barna libraries are released under the BSD-3 license

Copyright (c) 2010, Micha Sammeth
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The names of its contributors may be not used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

LIBRARIES
---------------
The 3rd party libraries used in the Flux Capacitor are released under the
following licenses:

License: 'GNU Lesser General Public License'
URL: http://www.gnu.org/licenses/lgpl.txt
  jsap 2.1
  jcommon 1.0.16
  jfreechart 1.0.13

License: 'LGPL 2.1'
URL: http://www.gnu.org/licenses/lgpl-2.1.html
  javassist 3.12.1.GA


License: 'The Apache Software License, Version 2.0'
URL: http://www.apache.org/licenses/LICENSE-2.0.txt
  gson 2.1
  guava r08
  commons-cli 1.2
  commons-math 2.2
  groovy-all 1.8.4
  xml-apis 1.0.b2

License: 'Mozilla Public License'
URL: http://www.mozilla.org/MPL/MPL-1.1.html
  itext 2.0.7

License: 'BSD-3'
URL: http://xstream.codehaus.org/license.html
  xstream 1.2.2

License: 'BSD'
URL: http://dom4j.sourceforge.net/dom4j-1.6.1/license.html
  dom4j 1.6.1

License: 'WTFL'
URL: http://sam.zoy.org/wtfpl/COPYING
  reflections 0.9.5

License: 'MIT'
URL: http://www.opensource.org/licenses/MIT
  slf4j-api 1.6.1
  slf4j-nop 1.6.1


License: 'Indiana University Extreme! Lab Software License'
URL: http://www.extreme.indiana.edu/dist/java-repository/xpp3/distributions/
  xpp3_min 1.1.3.4.O

The FluxDev-Team

Thasso and Micha
