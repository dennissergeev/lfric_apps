import re
import sys
from metomi.rose.upgrade import MacroUpgrade


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


class vn20_t85(MacroUpgrade):
    """Upgrade macro for ticket #85 by Chris Smith."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t85"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add new geostrophic_forcing namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:formulation",
            r"namelist:formulation" + "\n" + " (namelist:geostrophic_forcing)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add geostrophic_forcing setting to external_forcing namelist"""
        self.add_setting(
            config,
            ["namelist:external_forcing", "geostrophic_forcing"],
            ".false.",
        )
        """Add default data for geostrophic_forcing namelist"""
        self.add_setting(config, ["namelist:geostrophic_forcing"])
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "coordinate"], "'height'"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "heights"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_u"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_v"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "times"], "0.0"
        )

        return config, self.reports


class vn20_t358(MacroUpgrade):
    """Upgrade macro for ticket #358 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0_t85"
    AFTER_TAG = "vn2.0_t358"

    def upgrade(self, config, meta_config=None):
        # Commands From: components/driver/rose-meta/lfric-driver
        """
        Add element_order_h and element_order_v to namelist finite_element,
        replacing element_order
        """
        self.add_setting(
            config, ["namelist:finite_element", "element_order_h"], "0"
        )
        self.add_setting(
            config, ["namelist:finite_element", "element_order_v"], "0"
        )
        self.remove_setting(
            config, ["namelist:finite_element", "element_order"]
        )

        return config, self.reports


class vn20_t82(MacroUpgrade):
    """Upgrade macro for ticket #82 by Chris Smith."""

    BEFORE_TAG = "vn2.0_t358"
    AFTER_TAG = "vn2.0_t82"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add vapour_forcing namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = source + "\n" + " (namelist:vapour_forcing)"
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add vapour_forcing setting to external_forcing namelist"""
        self.add_setting(
            config, ["namelist:external_forcing", "vapour_forcing"], "'none'"
        )
        """Data for vapour_forcing namelist"""
        self.add_setting(config, ["namelist:vapour_forcing"])
        self.add_setting(
            config, ["namelist:vapour_forcing", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "times"], "0.0")

        return config, self.reports
