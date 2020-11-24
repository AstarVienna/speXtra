"""
Configuration for speXtra
"""

from astropy.config import ConfigItem, ConfigNamespace
import yaml

__all__ = ["Conf"]


class Config:
    def __init__(self, config_file):
        self.config_file = config_file
        with open(self.config_file) as f:
            self.meta = yaml.safe_load(f)
        for key in self.meta:
            setattr(self, key, self.meta[key])

    def set_param(self, name, value):
        d = {name: value}
        self.meta.update(d)
        self.__dict__.update(d)
        with open(self.config_file, "w") as f:
            yaml.dump(self.meta, f)

    def __repr__(self):
        return yaml.dump(self.meta, indent=4, sort_keys=False, default_flow_style=False)



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
