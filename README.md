# DYFF-1.0eh

EasyHybrid's extended DYFF force field for [pDynamo3](https://github.com/pdynamo/pDynamo3)
(this package needs [pdynamo3-eh](https://github.com/ferbachega/pdynamo3-eh), EasyHybrid's
pDynamo3 fork): an extension of pDynamo3's own generic, rule-based `dyff-1.0` with

- 175 additional atom types and 246 patterns (transition metals: a coordinate-geometry
  perceiver, η-bonded and high-spin offsets; aromatic P/Se/B, borane clusters),
- a biomolecular residue library (force-matched to DFT-level relaxed scans) for amino
  acids, nucleotides and common cofactors,
- CM5 partial charges (GFN1-xTB) and a small Lennard-Jones type for polar hydrogens,
- per-bond torsion refinement against relaxed GFN2-xTB scans,
- Hessian-derived (Seminario-style) bond/angle terms for metal sites.

**Status (2026-11):** just split out of
[EasyHybrid3](https://github.com/ferbachega/EasyHybrid3)'s `src/pdynamo/`. It does not need
the EasyHybrid GUI — only `pdynamo3-eh` and this package.

## Why "pdynamo" and not "dyffeh" as the package name

EasyHybrid itself (and every GUI window, the terminal, the automation API, ...) imports this
code as `pdynamo.dyff_edits`, `pdynamo.residue_library`, etc. Renaming it would have meant
touching every one of those call sites for no functional benefit, so this package keeps the
same name and relies on [PEP 420](https://peps.python.org/pep-0420/) implicit namespace
packages: EasyHybrid's own `src/pdynamo/` (the GUI-coupled bridge code that stays there,
never part of this package) and this package's `pdynamo/` merge into one `pdynamo` namespace
as long as both are on `sys.path` and neither has an `__init__.py` (neither does). Installing
this package standalone, with no EasyHybrid around, works exactly the same way — `pdynamo`
just resolves to this package's `pdynamo/` alone.

## Installing

```bash
git clone https://github.com/ferbachega/pdynamo3-eh.git
cd pdynamo3-eh/installation && python3 Install.py -f && cd -
git clone https://github.com/ferbachega/dyff-1.0eh.git
export PDYNAMO3_EH_HOME=$(pwd)/pdynamo3-eh      # or put it next to this repo as a sibling
export PYTHONPATH=$PDYNAMO3_EH_HOME:$(pwd)/dyff-1.0eh:$PYTHONPATH
```

xTB (for CM5 charges and torsion scans) and Open Babel (optional, native ligand perception
is the default — see `pdynamo/ligand_perception.py`) are external binaries, found on `PATH`
or at `$XTBHOME/bin/xtb` / `~/programs/xtb-*/bin/xtb`.

## Using it

```python
from pdynamo.dyff_edits import dyff_model, EH_DYFF_SET
from pdynamo import dyff_prep, residue_library

model = dyff_model(parameter_set=EH_DYFF_SET)   # dyff-1.0-eh
system.DefineMMModel(model)
```

See `examples/` in the EasyHybrid3 repository (`examples/dyff_prep_bachega/prepare_dyff_eh.py`)
for a full worked script: structure file in, DYFF-1.0eh system out, no GUI.

## Package layout

- `pdynamo/dyff_edits.py` — the DYFF-1.0eh parameter set switch (`dyff_model()`), the
  override mechanism that layers calibrated values on top of pDynamo3's own rule-based DYFF
  factories (see that module's own header for why: a pDynamo3 bug in reading explicit
  parameter-file rows), polar-H/TIP3P Lennard-Jones, metal bond/angle offsets.
- `pdynamo/dyff_prep.py` — whole-system DYFF preparation (metal site detection, side-chain
  completion, protonation).
- `pdynamo/residue_library.py` — the force-matched residue library: matching a system's
  residues to templates, charge assignment (CM5/RESP), per-bond torsion scanning.
- `pdynamo/metal_geometry.py` / `metal_refine.py` — transition-metal coordinate-geometry
  perception (pDynamo3's own typer only counts neighbours) and per-system Hessian refinement.
- `pdynamo/seminario.py` — Hessian-compliance (modified Seminario) bond/angle/out-of-plane
  fitting.
- `pdynamo/bond_library.py` / `dyff_bond_fit.py` / `torsion_fit.py` — the offline
  calibration pipeline (environment-code bond/angle/torsion libraries, per-bond fitting
  against relaxed xTB scans) that produced `pdynamo/data/*.json`, and the runtime fitting
  code it shares with the online (residue-library) path.
- `pdynamo/esp_charges.py` — CM5/RESP/ESP charges (GFN1/GFN2-xTB, ORCA).
- `pdynamo/env_codes.py` — the RDKit-free atom-environment code used to key the bond/angle
  libraries (a pure-Python port of the RDKit canonical-ranking algorithm the libraries were
  originally built with; `EASYHYBRID_ENV_CODES=rdkit` restores the RDKit path, `pip install
  dyff-1.0eh[rdkit]`).
- `pdynamo/ligand_perception.py` — native (RDKit/Open-Babel-free) ligand chemistry
  perception: hydrogen placement, Wang & Case bond orders, pH-dependent protonation.
- `pdynamo/bond_order_perception.py` — the Wang & Case bond-order-perception algorithm
  itself (vendored 2026-11 from EasyHybrid's `graphics_engine`/vismol: pure graph algorithm,
  stdlib only, no graphics dependency at all — `ligand_perception.py` and
  `ligand_chemistry.py` used to reach it through a GUI-submodule-relative path trick that
  only worked when the whole graphics engine happened to be present, and otherwise failed
  *silently* into a cruder fallback; found and fixed while testing this package standalone).
- `pdynamo/ligand_chemistry.py` / `structure_prep.py` — force-field-agnostic ligand
  chemistry and whole-system analysis/building (split out of EasyHybrid's `pdynamo/opls/`
  2026-11: that module mixed this generic part with OPLS-AA-specific atom typing; DYFF's
  preparation only ever needed the generic part).
- `pdynamo/mol2_import.py` — fallback MOL2 reader for charged aromatics pDynamo3's own
  reader rejects.
- `pdynamo/subsystem.py` — capped subsystems (torsion-scan fragments, selection refinement).
- `pdynamo/param_logs.py` — keeps the original xtb/ORCA outputs of a parametrization.
- `pdynamo/data/` — the bond/angle/hessian/residue libraries (JSON) and the `dyff-1.0-eh`
  atom-type/pattern set itself (`data/forcefields/dyff/dyff-1.0-eh/`, in pDynamo3's own
  native parameter-set format).
- `pdynamo/pdb_components_extra/` — PDB component library entries pDynamo3's own library
  lacks (CYM, LYN protonation variants).

## Tests

```bash
pip install -e .[dev]
pytest tests
```
