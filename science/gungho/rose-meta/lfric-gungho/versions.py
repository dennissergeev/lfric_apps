import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version20_21 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn21_t83(MacroUpgrade):
    """Upgrade macro for ticket #83 by Chris Smith."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t83"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add theta_relax namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"\(namelist:temp_tend_data\)",
            r"(namelist:temp_tend_data)" + "\n" + " (namelist:theta_relax)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add theta_relaxation setting to external_forcing namelist"""
        self.add_setting(
            config, ["namelist:external_forcing", "theta_relaxation"], ".false."
        )
        """Data for theta_relax namelist"""
        self.add_setting(config, ["namelist:theta_relax"])
        self.add_setting(
            config, ["namelist:theta_relax", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:theta_relax", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:theta_relax", "number_heights"], "1"
        )
        self.add_setting(config, ["namelist:theta_relax", "number_times"], "1")
        self.add_setting(
            config, ["namelist:theta_relax", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:theta_relax", "times"], "0.0")
        self.add_setting(config, ["namelist:theta_relax", "timescale"], "1.0")
        return config, self.reports


class vn21_t476(MacroUpgrade):
    """Upgrade macro for ticket #476 by Dave Case."""

    BEFORE_TAG = "vn2.1_t83"
    AFTER_TAG = "vn2.1_t476"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Set configure_segments in physics namelist to be false to ignore segmentation"""
        self.add_setting(config, ["namelist:physics", "bl_segment"], "0")
        self.add_setting(config, ["namelist:physics", "gw_segment"], "0")
        self.add_setting(config, ["namelist:physics", "ussp_segment"], "0")

        return config, self.reports


class vn21_t663(MacroUpgrade):
    """Upgrade macro for ticket #663 by Ian Boutle."""

    BEFORE_TAG = "vn2.1_t476"
    AFTER_TAG = "vn2.1_t663"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Commands From: rose-meta/um-aerosol
        self.add_setting(config, ["namelist:aerosol", "murk_lbc"], ".false.")

        return config, self.reports
