"""
Configuration for speXtra
"""

from astropy.config import ConfigItem, ConfigNamespace

__all__ = ["Conf"]


class Conf(ConfigNamespace):
    """Configuration parameters for spextra.
    TODO: Write default configuration to an archive so it can be _really_ changed in runtime
    """

    data_dir = ConfigItem(None,
                          "Directory where spextra will store and read downloaded data "
                          "resources. "
                          "If None, ASTROPY_CACHE_DIR/spextra is created and used. "
                          "Example of user defined: data_dir = '/home/user/data/spextra' ",
                          cfgtype='string(default=None)')

    database_url = ConfigItem("https://homepage.univie.ac.at/miguel.verdugo/database/",
                              "URL address of the database",
                              cfgtype='string(default=None)')

    remote_timeout = ConfigItem(30.0, "Remote timeout in seconds.")
