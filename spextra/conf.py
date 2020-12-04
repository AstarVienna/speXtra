"""
Configuration for speXtra
"""

from astropy.config import ConfigItem, ConfigNamespace
import yaml

__all__ = ["Conf"]


class Config:
    def __init__(self, data_dir=None, database_url=None, remote_timeout=None):

        self.config_file = __config_file__

        with open(self.config_file) as f:
            self.meta = yaml.safe_load(f)
            for key in self.meta:
                setattr(self, key, self.meta[key])

        if data_dir is not None:
            self.set_param("data_dir", data_dir)

        if database_url is not None:
            self.set_param("database_url", database_url)

        if remote_timeout is not None:
            self.set_param("remote_timeout", remote_timeout)

    def set_param(self, name, value):
        d = {name: value}
        self.meta.update(d)
        self.__dict__.update(d)
        with open(self.config_file, "w") as f:
            yaml.dump(self.meta, f, sort_keys=False)

    def __repr__(self):
        return yaml.dump(self.meta, indent=4, sort_keys=False, default_flow_style=False)

    def get_rootdir(self):
        # use the environment variable if set
        data_dir = os.environ.get('SPEXTRA_DATA_DIR')

        # otherwise, use config file value if set.

        if data_dir is None:
            data_dir = self.data_dir

        # if still None, use astropy cache dir (and create if necessary!)
        if data_dir is None:
            data_dir = os.path.join(get_cache_dir(), "spextra")
            if not os.path.isdir(data_dir):
                if os.path.exists(data_dir):
                    raise RuntimeError("{0} not a directory".format(data_dir))
                os.mkdir(data_dir)

        return data_dir

    def get_database_url(self):

        """
        Check if the database is reachable

        Returns
        -------
        the database_location
        """
        loc = self.database_url
        try:
            assert is_url(loc)
        except AssertionError as error:
            print(error)
            #  warnings.warn("Database address not reachable", loc)
            pass

        return loc


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
