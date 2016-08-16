ASTALAVISTA
===========

DESCRIPTION
-----------

Astalavista is a program to analyze splicing structures of gene annotations.
In a standard application it retrieves a set of so-called events according to
a consistent definition.

REQUIREMENTS
------------

* java 1.6 or higher installed

* not more than 1Gb of RAM

* a gene annotation that is to be analyzed

CHANGES
---------------

AStalavista 4.0
    * released astafunk tool
    * [BARNA-386] - astalavista scorer always complains of pre-existing 'h' flag

AStalavista 3.2
    * [BARNA-304] - Astalavista starter gives an error on startup
    * [BARNA-305] - Scorer specific parameters --vcf and --sp are not recognized
    * [BARNA-311] - AStalavista always outputs CDS events

AStalavista 3.1
    * [BARNA-133] - Connect AStalavista tool
    * [BARNA-227] - Implement Site Scoring
    * [BARNA-270] - AStalavista option to output ALL variations of the exon-intron structure

GETTING STARTED
---------------

* unpack

 If you can read this, you obviously managed to unpack Astalavista :) Welcome

 You will find the executable in the bin/ folder

    astalavista.bat     Windows

    astalavista         UNIX clones and Mac OSX


* Run

    You can take a look at the command line parameters and available tools using

     astalavista --help

    The Flux Capacitor is one of the tools in the Flux Package, further information
    on each of the tools can be retrieved by the command pattern

     astalavista --t <toolname> --help


* Parameters

    Astalavista parameters can be specified via the command line or in a separate
    parameter file with the -p command line parameter (see next section)

    To get a list and descriptions for all available parameters, use

    astalavista --printParameters

    to create a file with an exhaustive list of parameters and their default values,
    pipe the output to a file

    astalavista --printParameters > myparameters.par

* Memory

    In case you run into out of memory issues you can increase the memory size used
    by Astalavvista using the environment variable FLUX_MEM, for example:

    # this sets the upper limit to 2 gig
    export FLUX_MEM="2G"; astalavista -p ....


    # this also sets the upper limit to 2 gig
    export FLUX_MEM="2048M"; astalavista -p ....

* The program homepage

	http://astalavista.sammeth.net

* Especially, have a look at the documentation of AStalavista at

	http://sammeth.net/confluence/display/ASTA/Home

* If you encounter any bugs, use the bugtracking system at

	hhttp://sammeth.net/jira

  First check for any known bugs, if you don't find your issue described, log
  in and create a new issue.

* If you have any questions, discuss them via the forum

	http://sammeth.net/confluence/display/ASTA/Appendix+D+-+Forum

or leave me an email.

LICENSE
---------------
The Barna libraries containing the Astalavista program are released under the BSD-3 license

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
The 3rd party libraries used of Astalavista are released under the
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

Micha (on behalf of all Astalavista developers)
