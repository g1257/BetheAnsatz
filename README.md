# Quick Start
 
# Disclaimer and Licensing
 
BetheAnsatz is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
BetheAnsatz is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with BetheAnsatz. If not, see <http://www.gnu.org/licenses/>.
The full software license for BetheAnsatz version 1.0.0 
can be found in
file COPYING. 

# Please cite this work

BetheAnsatz is a free and open source Bethe ansatz code.
The full software license for BetheAnsatz version 0.1
can be found in
file COPYING. 
You are welcomed to use it and publish data 
obtained with BetheAnsatz. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

# Hash of the latest commit 

Hash of the latest commit is also posted at
https://web.ornl.gov/~gz1/hashes.html

# Building and Running BetheAnsatz

## Required Software

* GNU C++
* The LAPACK and BLAS libraries
* The GSL library
* PsimagLite (see below)

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script) 

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
the LAPACK and BLAS libraries, the gsl library, make, perl, and git 
if you don't have them.

2. Issue

    cd someDirectory/

    git clone https://github.com/g1257/PsimagLite.git

    git clone https://github.com/g1257/BetheAnsatz.git

3. Compile PsimagLite

    cd PsimagLite/lib/

    make -f Makefile.sample

    cd ../../

4. Now issue

    cd BetheAnsatz/src

    make

5. You can run it with TBW.
 
