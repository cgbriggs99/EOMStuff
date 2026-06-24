#
# @BEGIN LICENSE
#
# general_cc by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util


def run_general_cc(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    general_cc can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('general_cc')

    """
    psi4.core.reopen_outfile()
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes hereN
    # psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)
    
    molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    
    if molecule.schoenflies_symbol() != 'c1':
        psi4.core.print_out("This code must be run in C1 symmetry! Changing symmetry.");
        molecule = molecule.clone()
        molecule.reset_point_group('c1')
        molecule.update_geometry()
        
    psi4.core.print_out("Creating basis sets...\n");
    
    scf_aux_basis = psi4.core.BasisSet.build(molecule, key="DF_BASIS_SCF", target=psi4.core.get_option("SCF", "DF_BASIS_SCF"),
                                         fitrole="JKFIT", other=psi4.core.get_global_option("BASIS"))
    
    mp2_aux_basis = psi4.core.BasisSet.build(molecule, key="DF_BASIS_MP2", target=psi4.core.get_option("DFMP2", "DF_BASIS_MP2"),
                                         fitrole="RIFIT", other=psi4.core.get_global_option("BASIS"))
    
    psi4.core.timer_on("My RHF code")
    
    ref_wfn = psi4.core.Wavefunction.build(molecule, psi4.core.get_global_option("BASIS"))
    
    psi4.core.print_out("Setting basis sets...");
    
    if psi4.core.get_global_option("SCF_TYPE").upper() == "DF":
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)
        
    if psi4.core.get_global_option("MP2_TYPE").upper() == "DF":
        ref_wfn.set_basisset("DF_BASIS_MP2", mp2_aux_basis)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    cgbriggs_wfn = psi4.core.plugin('general_cc.so', ref_wfn)
    
    psi4.set_variable("CURRENT ENERGY", cgbriggs_wfn.energy())
    
    psi4.core.timer_off("My RHF code")

    if kwargs.get('return_wfn', False):
        return cgbriggs_wfn.energy(), cgbriggs_wfn
    return cgbriggs_wfn.energy()


# Integration with driver routines
psi4.driver.procedures['energy']['general_cc'] = run_general_cc


def exampleFN():
    # Your Python code goes here
    pass
