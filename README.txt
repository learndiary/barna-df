1. What is the Barna-df project?

"Barna-df" project is a dirty-fixed of the original Barna project on
"https://confluence.sammeth.net/".

There are "Flux Simulator", "Astalavista" and "Flux Capacitor" these three
softwares in Barna project.
 
It seems the original Barna project has broken. But, there are still people use
it. so I uploaded a dirty fix of this project as "barna-df (barna dirty fix)"
project. There is no any change of original working code, only fixed the
invalid dependent package repositories or packages to make the source building
can be done.

Everything else such as documentation, issue tracker and etc. please see
original Barna project home on "https://confluence.sammeth.net/".

2. How to build project from source on Linux?

# Prepare environment
Please make sure that you have a full jdk version 1.7 or 1.8 installed, i.e.,
not just a runtime environment without the 'javac' compiler binary. You can
check the presence of the compiler commnd in JAVA_HOME/bin. Also make sure that
the JAVA_HOME environment variable is set to the root folder of your Java
installation, by

export JAVA_HOME=/path/to/java

# Get the source
$> git clone https://github.com/learndiary/barna-df.git

# Applying the dirty-fixed patch
$> cd barna-df
$> patch -p1 < patch/barna-df.patch

# Building Flux Simulator
$> cd barna.simulator
$> ../gradlew dist

When executed the first time, gradle will download a lot of dependencies, get a
hot beverage, this might take a moment. The gradle script will compile the code
and all dependencies, and finally a .tag.gz and a .zip file are created. You
can find them in the build/distributions/ folder.These are the created bundles
that contain the Flux Simulator build from source.

Astalavista and Flux Capacitor below creates the similiar built binary result
as Flux Simulator above.

# Building Astalavista
$> cd ../barna.astalavista
$> ../gradlew dist

# Building Flux Capacitor
$> cd ../barna.capacitor
$> ../gradlew dist

3. How to use it?

For example, we need use Flux Simulator built above.

$> cd ../barna.simulator/build/distributions
$> tar -xf flux-simulator-1.2.3.tgz
$> cd flux-simulator-1.2.3/bin
$> ./flux-simulator --version # This will show the version information as below

[INFO] Flux-Simulator v1.2.3 (Flux Library: 1.30)

Flux-Simulator
Version 1.2.3
Barna Library 1.30
-----------------------------------------------
Build Date Fri Aug 19 10:23:32 CST 2022
Build Version 2cb4c7acd80c4347d1e03f10a68110c1140867b1
Build Branch

It is the similiar way to use AStalavista and Flux Capacitor as Flux Simulator
above.
 
Below is the content of README.txt file of original Barna project.

######################################################

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
