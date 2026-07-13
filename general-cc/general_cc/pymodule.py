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


def generate_set(level, l, n, k, prev = None) :
    out = []
    bound1 = min(level, l)
    if prev is not None :
        bound1 = min(bound1, prev[0])
    for m1 in range(0, bound1 + 1) :
        bound2 = min(m1, n)
        if prev is not None and m1 == prev[0]:
            bound2 = min(bound2, prev[1])
        for m2 in range(0, bound2 + 1) :
            bound3 = min(m2, k)
            if prev is not None and m1 == prev[0] and m2 == prev[1] :
                bound3 = min(bound3, prev[2])
            for m3 in range(0, bound3 + 1) :
                add = [(m1, m2, m3)]
                if l - m1 == 0 and n - m2 == 0 and k - m3 == 0 :
                    out.append(add)
                else :
                    new_list = generate_set(level, l - m1, n - m2, k - m3)
                    for elem in new_list :
                        out.append(add + elem)
    return out



def generate_density(level, particles) :
    out = {}
    for num_part_an in range(0, particles + 1) :
        for num_hole_an in range(0, particles + 1) :
            elements = {}
            n = num_part_an + num_hole_an
            k = num_part_an
            num_part_cr = particles - num_hole_an
            num_hole_cr = particles - num_part_an
            m = min(num_part_cr, num_hole_cr)
            for lam in range(max(1, m + n - 2 * min(n - k, k)), level + 1) :
                lam_excess = lam - m - n + 2 * min(n - k, k)
                
                elements[(lam, lam_excess)] = generate_set(level, lam_excess, n, k)
            out[f"{num_part_cr}pc,{num_hole_an}ha,{num_hole_cr}hc,{num_part_an}pa"] = elements
    return out



