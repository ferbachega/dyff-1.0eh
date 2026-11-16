#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple_prepare.py -- the fewest steps to get a DYFF-1.0eh system.
#
#  For the fully-commented, step-by-step version (reports, metal sites,
#  minimization, CLI options), see examples/dyff_prep_bachega/prepare_dyff_eh.py.
#  This one is deliberately short: a structure file in, a DYFF-1.0eh system out.
#
#  Usage:
#      python3 simple_prepare.py structure.pdb
#      python3 simple_prepare.py                 # the bundled 1BX4.pdb
#

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PDYNAMO3_EH = os.environ.get("PDYNAMO3_EH_HOME") or os.path.join(os.path.dirname(REPO), "pdynamo3-eh")
for path in (REPO, PDYNAMO3_EH):
    if path not in sys.path:
        sys.path.insert(0, path)

from pBabel import ImportSystem                          # noqa: E402
from pCore import Pickle                                 # noqa: E402
from pMolecule.NBModel import NBModelCutOff               # noqa: E402

from pdynamo import dyff_prep as DP                       # noqa: E402
from pdynamo import residue_library as RL                 # noqa: E402
from pdynamo.dyff_edits import dyff_model, soften_polar_hydrogens, water_lj , \
                               apply_guanidinium_correction                 # noqa: E402


def prepare_dyff(path, ph=7.4):
    """ Structure file -> pDynamo System with the DYFF-1.0eh MM model, charges
        and a non-bonded model defined. """

    # 1. Read the file and work out its chemistry: protonation, hydrogens,
    #    metal sites (System > Prepare DYFF System in the GUI).
    raw = ImportSystem(path, log=None)
    info = DP.analyze(raw)
    DP.apply_site_protonation(info)
    system, prep = DP.build(info, complete_truncated=True, ph=ph)

    # 2. Residue library: charges and bonded terms for every residue it
    #    recognizes (amino acids, nucleotides, water, ions, ...).
    plan = RL.plan(system)
    charges = dict(plan["charges"])
    edits = plan["edits"]

    # 3. Anything left (ligands, unusual residues) gets CM5 charges from
    #    GFN2-xTB. hessian=False/torsions=False: charges only, the fastest
    #    option (see prepare_dyff_eh.py for the full xTB Hessian/torsion fit).
    if plan["rest_atoms"]:
        rest = RL.rest_parameters(system, plan["rest_atoms"], hessian=False, torsions=False)
        charges.update(rest["charges"])

    # 4. Build the DYFF-1.0eh MM model and set the charges.
    apply_guanidinium_correction(system)
    system.DefineMMModel(dyff_model("dyff-1.0", edits, prep.get("type_overrides")), log=None)
    for i, q in charges.items():
        system.mmState.charges[i] = q
    # . not normalized to an exact integer total here (prepare_dyff_eh.py spreads the
    #   rounding residual over every atom) -- usually off by < 1e-3 e, fine to ignore
    #   for a quick look, but do that normalization for production work.
    soften_polar_hydrogens(system)      # small LJ type for polar H (needs charges)
    water_lj(system)                    # TIP3P oxygen LJ for water

    # 5. A non-bonded model, so the system has an energy.
    system.DefineNBModel(NBModelCutOff.WithDefaults())
    return system


if __name__ == "__main__":
    structure = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "dyff_prep_bachega", "1BX4.pdb")
    system = prepare_dyff(structure)
    print("{} atoms, E = {:.2f} kJ/mol".format(len(system.atoms), system.Energy(log=None)))
    out = os.path.splitext(os.path.basename(structure))[0] + "_dyff_eh.pkl"
    Pickle(out, system)
    print("saved:", out)
