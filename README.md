                       +------------------------------+
                       |     The Elk FP-LAPW Code     |
                       +------------------------------+

This code is distributed under the terms of the GNU General Public License.
See the file COPYING for license details.

    Elk is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Elk is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with Elk. If not, see http://www.gnu.org/licenses/. 

Elk can be compiled by first running "./setup" in this directory followed by
"make all". This will compile the main code as well as several auxiliary
programs. For optimal performance we strongly recommend that you tune the
Fortran compiler options in the file "make.inc" and use machine-optimised
BLAS/LAPACK libraries. Setting the OpenMP options of your compiler will enable
Elk to run in parallel mode on multiprocessor systems.

A test suite is available: entering "make test" will check the output of your
executable against a standard set. This may take some time to complete.

Auxiliary programs include "spacegroup" for producing crystal geometries from
spacegroup data, and "eos" for fitting equations of state to energy-volume data.

Elk is updated regularly with new features and bug fixes. Features not listed as
"experimental" may be used for production but, as with any code, please check
the consistency of your results carefully.

--------------------------------------------------------------------------------
J. K. Dewhurst, S. Sharma
L. Nordstrom,  F. Cricchio, F. Bultmark
E. K. U. Gross

Notes on the GW module
======================

The current GW module includes two sub-modules that allow one to compute the 
quasi-particle energies in (1) real-frequency domain, and (2) Matsubara-time domain.
The details of the implementations can be found at::

1. I.-H. Chu, A. Kozhevnikov, T. C. Schulthess, and H.-P. Cheng, "All-electron 
GW quasiparticle band structures of group 14 nitride compounds", J Chem. Phys. 
141, 044709 (2014), doi:http://dx.doi.org/10.1063/1.4890325.

2. I.-H. Chu, J. P. Trinastic, Y.-P. Wang, A. G. Eguiluz, A. Kozhevnikov, T. C. 
Schulthess, and H.-P. Cheng, "All-electron self-consistent GW in the Matsubara-time 
domain: Implementation and benchmarks of semiconductors and insulators", Phys. Rev. B 
93, 125210 (2016), doi:https://doi.org/10.1103/PhysRevB.93.125210
	
