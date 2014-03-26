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

FluxCapacitor 1.6
    * [BARNA-328] - When writing to stdout, the flux prints its version information which ends up in the GTF output
    * [BARNA-366] - Joint deconvolution of paired-end reads
    * [BARNA-367] - Simplification of delta values

FluxCapacitor 1.5.2
    * [BARNA-340] - implement parameters for annotation sources to be included or excluded in profiling
    * [BARNA-341] - stabilize profiling in the presence of highly abundant transcripts added to the transcriptome annotation

FluxCapacitor 1.5.1
    * [BARNA-339] - Array overflow in profiling

FluxCapacitor 1.5
    * [BARNA-334] - Replace Attributes with Mapping properties
    * [BARNA-336] - SAMException - Value for tag XT is not a Character
    * [BARNA-337] - Improve profiles stability
    * [BARNA-338] - Collapse transcripts that are exactly identical in terms of genomic coordinates of their exons

FluxCapacitor 1.4
    * [BARNA-329] - AnnotationMapper maps genomic mappings spanning all-intronic edges
    * [BARNA-330] - Fraction Inconsistency

FluxCapacitor 1.3
    * [BARNA-28] - Exit with error 1 in case of error
    * [BARNA-273] - Implement stranded annotation mapping for BAM input files
    * [BARNA-296] - OutOfMemoryError for some loci while getting the sorted iterator for BAM files
    * [BARNA-301] - NullPointerException while iterating loci in BAM files
    * [BARNA-308] - Capacitor does not start with single end data
    * [BARNA-309] - Capacitor creates lots of threads when running on BAM files
    * [BARNA-314] - flux capacitor error:java.io.IOException: Pipe closed
    * [BARNA-321] - NullPointer when trying to set LP output
    * [BARNA-323] - Add SAM validation stringency parameter to SAMMappingSortediterator
    * [BARNA-315] - Parameters for Mapping Filtering
    * [BARNA-325] - Make sort-in-ram default
    * [BARNA-326] - Remove the no longer needed NR_READS_MAPPED parameter
    * [BARNA-327] - Catch ParameterException in main method to hide stack trace

FluxCapacitor 1.2.4
    * [BARNA-289] - Change NO_DECOMPOSE parameter to DECONVOLUTE
    * [BARNA-297] - Stats file lock might fail and then no stats are written
    * [BARNA-298] - Capacitor prints NaN for read count and RPKM
    * [BARNA-310] - incorrect junction and intron coordinates

FluxCapacitor 1.2.3
    * [BARNA-253] - Add Differential-Expression Tool to the Capacitor
    * [BARNA-264] - Differential Expression result location misses the : between start and end
    * [BARNA-265] - error index out of bounds fluxcapacitator
    * [BARNA-266] - Capacitor crashed with BED input file
    * [BARNA-267] - The stats file writer tries to lock the stats file, but that might fail
    * [BARNA-268] - Capacitor throw an NPE for Ensmeble annotation
    * [BARNA-269] - Capacitor throw NPE when trying to ouput coverage stats

FluxCapacitor 1.2.2
    * [BARNA-263] - BAM ID parsing error results in NPE

FluxCapacitor 1.2.1
    * [BARNA-261] - Capacitor complains about read descriptor for bam files

FluxCapacitor 1.2
    * [BARNA-215] - Sequencing fails with ArrayIndexOutOfBoundsException
    * [BARNA-218] - Replace BEDObject2 with Mapping
    * [BARNA-235] - Flux-simulator does not accept SIZE_DISTRIBUTION N(300,30)
    * [BARNA-236] - Single-end stranded data give 0 RPKM for trancripts with few reads supporting them
    * [BARNA-239] - Parameter Schema relative path parser extends Windows paths
    * [BARNA-240] - Capacitor runs out of Heap Space with some BAM files
    * [BARNA-241] - net.sf.samtools.SAMFormatException: SAM validation error
    * [BARNA-246] - UnixSorter line interceptor called for each merge step instead of just the final one
    * [BARNA-247] - Capacitor still wants the -p <parameterfile> option even if I specify all mandatory paramters on the command line
    * [BARNA-248] - Capacitor complains about READ_DESCRIPTOR when given aunindexed BAM file
    * [BARNA-251] - Parameter file parser does not undestand # comments after the parameter
    * [BARNA-256] - When an error occurs solving the linear system, the matrix is written on stderr and an exception is thrown even if the deconvolution continuesi
    * [BARNA-237] - Non deterministic number of reads and RPKM on Capacitor old versions
    * [BARNA-238] - Wrong number of reads detected in some BAM input files
    * [BARNA-242] - When counting Splice Junctions and Introns, output the correct GeneID of a locus if multiple Genes are contained
    * [BARNA-245] - Scanning mapping file with BAM input could take a very long time
    * [BARNA-252] - Flux launcher fails on win32 when trying to allocate more than 1.5G for JVM heap
    * [BARNA-257] - Error in read ID: could not parse read identifier in indexed bam generated by GEMtools


FluxCapacitor 1.1
    * [BARNA-13] - UnixStreamSorter shouldn't close streams
    * [BARNA-22] - Synchronization between Log.setInteractive() and CommandLine.confirm()
    * [BARNA-25] - LineComparator Caching
    * [BARNA-114] - Capacitor stats are not written properly sometimes
    * [BARNA-119] - Line length is used instead of the mapping length for retrieving coverage stats within the learn() method
    * [BARNA-131] - Capacitor Tests run into endless runs when GTF reader is not properly initialized (NULL pointer)
    * [BARNA-190] - FC need full path for output file on command line
    * [BARNA-194] - Reading gene sequences close to the chromosome border
    * [BARNA-199] - FC output wrong strand information for countings
    * [BARNA-200] - The gradle builds fail on windows because we call git directly
    * [BARNA-207] - FluxCapacitor deletes existing uncompressed BED file
    * [BARNA-208] - FluxCapacitor does not use sorted file for unsorted input if it already exists
    * [BARNA-210] - Null pointer exception when calling flux -t {TOOL} with no parameters
    * [BARNA-211] - LPSolverLoader assumes that libraries in the destination folder are correct, if identical by name
    * [BARNA-18] - ParameterSchema with multiple Settings instances
    * [BARNA-68] - Refactor FluxSimulatorSettings.RelativePathParser to a more general class outside of the Simulator settings
    * [BARNA-87] - Check SyncIOhandler2 for removal
    * [BARNA-193] - Gene Identifier preserved in GTF wrapper
    * [BARNA-209] - Move KEEP_SORTED to File parameter
    * [BARNA-212] - Clean up class GraphLPsolver
    * [BARNA-214] - Unit test and integration tests do not delete all files
    * [BARNA-138] - Implement BAM Reader to read BAM files

FluxCapacitor 1.0.2
    - BARNA-188	Check Capacitor number of reads for transcript
    - BARNA-186	Update BED wrapper for CASAVA 1.8 read descriptor
    - BARNA-160	FC should accept multi-split mappings for annotation mapping
    - BARNA-149 Count reads to splice junctions
 	- BARNA-140	Too many open files
 	- BARNA-132 Count reads to all-intronic regions
    - BARNA-168 - ByteArrayCharSequence.complement() cannot handle IUPAC ambiguities

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

	http://sammeth.net/confluence/display/FLUX/Home

* If you encounter any bugs, use the bugtracking system at

	http://sammeth.net/jira
	
  First check for any known bugs, if you don't find your issue described, log
  in and create a new issue.

* If you have any questions, discuss them via the forum 
	
	http://sammeth.net/confluence/display/FLUX/Appendix+D+-+Forum
	
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
