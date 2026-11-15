"""
Shared pytest setup for the dyff-1.0eh test suite.

This package depends on pDynamo3 (pdynamo3-eh, EasyHybrid's fork -- see this
repository's README). Point PDYNAMO3_EH_HOME at your checkout of it; if unset,
a sibling ../pdynamo3-eh directory (next to this repository) is tried, which
is how the two are laid out side by side during development.
"""
import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_TESTS_DIR)

_PDYNAMO3_EH = os.environ.get("PDYNAMO3_EH_HOME") or os.path.join(
    os.path.dirname(_REPO_ROOT), "pdynamo3-eh"
)

for _p in (_REPO_ROOT, _PDYNAMO3_EH):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)
