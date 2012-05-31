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

FluxSimulator 1.0.1
    - BARNA-166 fix the error model loading

FluxSimulator 1.0
    - Fixed issue with Java 1.7
    - Fixed issue with the error model
    - Updated startup script to flux-simulator and set a default tool
    - Switched to 3G of RAM by default

FluxSimulator RC5
    -fixed Issue Issue 68: TMP_DIR parameter ignored - the user applied tmp folder should be used and
     debug tool added, just in case.

FluxSimulator RC4
    -fixed problem with disabling the PCR distriubution
    -fixed windows startup problems
    -fixed BED sorter issue

FluxSimulator RC3

    -fixed issue #58 - made sure all buffers are flushed before moving files
    -fixed issue #61 - temp directory is passed to FileUtils
    -fixed issue #60 - added documentation
    -fixed issue #48 - fastq output and new error model - you need an error model for for fasta/fastq now
                       NOTE : the FASTQ paramter is called FASTA now. If set to true, a fasta file will
                       be generated. If in addition ERR_FILE is specified, a fastq file is generated,
                       where the actual sequences incorporate errors.
    -added the ability to specify size distributions using normal distributions with mean and sd (i.e. N(200,20))
    -added PCR with a combined GC filtering step
    -added quality error models and a tool to create such models (currently only from GEM file)



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

    flux-simulator.bat     Windows

    flux-simulator         UNIX clones and Mac OS X


* Run

    You can take a look at the command line parameters and available tools using

     flux-simulator --help

    To start a specific tool, i.e. the error model generator, and see its parameters, use

     flux-simulator --t <toolname> --help


* Parameters

    The FluxSimulator reads its configuration from a parameter file specified with
    the -p command line parameter (see next section)

    To get a list and descriptions for all available parameters, use

    flux-simulator --printParameters

    to create a file with an exhaustive list of parameters and their default values,
    pipe the output to a file

     flux-simulator --printParameters > myparameters.par


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

    The Simulator documentation contains examples including demo data:

    http://barnaserver.com/confluence/display/SIM/4+-+Example+Simulations

* Memory

    In case you run into out of memory issues, you can increase the memory size used by
    the simulator using the environment variable FLUX_MEM, for example:

    # this sets the upper limit to 2 gig
    export FLUX_MEM="2G"; flux-simulator -p ....


    # this also sets the upper limit to 2 gig
    export FLUX_MEM="2048M"; flux-simulator -p ....

* All pages related to the program are reachable from the program homepage.
  Especially, have a look at the documentation of the FLUX SIMULATOR at

	http://barnaserver.com/confluence/display/SIM

* If you encounter any bugs, use the bugtracking system at

	http://barnaserver.com/jira/browse/BARNA
	
  First check for any known bugs, if you don't find your issue described, log
  in and create a new issue.


or leave me an email.

LICENSE
---------------
The Flux Simulator and all Barna libraries are released under the BSD-3 license

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
The 3rd party libraries used in the Flux Simulator are released under the
following licenses:

License: 'GNU Lesser General Public License'
URL: http://www.gnu.org/licenses/lgpl.txt
  jsap 2.1
  jcommon 1.0.0
  jfreechart 1.0.0


License: 'The Apache Software License, Version 2.0'
URL: http://www.apache.org/licenses/LICENSE-2.0.txt
  guava r08
  commons-cli 1.2
  jdbm 2.3
  commons-math 2.2
  groovy 1.8.4
  xml-apis 1.0.b2

License: 'Mozilla Public License 1.1'
URL: http://www.mozilla.org/MPL/MPL-1.1.html
  itext 2.0.7
  javassist 3.12.1.GA

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

License: 'Indiana University Extreme! Lab Software License'
URL: http://www.extreme.indiana.edu/dist/java-repository/xpp3/distributions/
    xpp3_min 1.1.3.4.O

The FluxDev-Team

Micha and Thasso
