import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn31_tXXXX(MacroUpgrade):
    """Upgrade macro for ticket #XXXX: random spherical harmonic internal
    flux forcing (equation 6 of Showman, Tan & Zhang, 2019, ApJ,
    arXiv:1807.08433)."""
    # TODO: replace tXXXX/vn3.1_tXXXX with the real ticket number.

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_tXXXX"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/jules_interface/rose-meta/jules-lsm
        nml = "namelist:specified_surface"
        self.add_setting(config, [nml, "internal_flux_forcing_amplitude"], "0.0")
        self.add_setting(config, [nml, "internal_flux_forcing_wavenumber"], "20")

        return config, self.reports
