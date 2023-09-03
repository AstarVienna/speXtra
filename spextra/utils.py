from os import devnull  # what is this?
import warnings
from contextlib import closing

# TODO: we should really use something better here...
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import socket
from pathlib import Path

import numpy as np
import astropy.units as u
# TODO: tqdm ffs
from astropy.utils.console import ProgressBarOrSpinner
from astropy.config import get_cache_dir
import yaml


__all__ = ["SpextraError",
           "Config",
           "download_file",
           "dict_generator",
           "download_svo_filter"]

ARRAY_LIKE = (np.ndarray, list, tuple)

__pkg_dir__ = Path(__file__).parent
__data_dir__ = __pkg_dir__ / "data"
__config_file__ = __data_dir__ / "config.yml"


class SpextraError(Exception):
    """Base class for errors in this Package."""


class Config:
    """
    Set up the configuration for the database.

    This class reads a yaml file in the directory data_dir that contains the
    default configuration of the database that contain default values for the
    following attributes:

    database_url: The location of the remote database
    data_dir: The location of the database in the local hard disk, where the
        files are saved. Default is Null (None) which then lead to
        .astropy/cache/spextra
    remote_timeout: The time before timeout

    Normally, the user will not need to call this class, except when chaning

    If the user needs to change these values the class can be simply called as

    >> Config(data_dir="path_to_new_data_directory")

    and the new path will be used from that moment on
    """

    def __init__(self, data_dir=None, database_url=None, remote_timeout=None):

        self.config_file = __config_file__

        with self.config_file.open(encoding="utf-8") as file:
            self.meta = yaml.safe_load(file)

        if self.meta["data_dir"] is not None:
            self.data_dir = Path(self.meta["data_dir"])
        else:
            self.data_dir = None
        self.database_url = self.meta["database_url"]
        self.remote_timeout = self.meta["remote_timeout"]
        self.default_data_dir = Path(get_cache_dir(),
                                     self.meta["default_data_dir"])

        if data_dir is not None:
            self._set_param("data_dir", Path(data_dir))

        if database_url is not None:
            self._set_param("database_url", database_url)

        if remote_timeout is not None:
            self._set_param("remote_timeout", remote_timeout)

    def get_data_dir(self) -> Path:
        """
        Retrieve the current data directory.
        
        Where the data from the database is stored.
        Check whether the directory exists and if not it will create it.
        Raises an error if not possible.

        Returns
        -------
        path: pathlib.Path
            Path to the local data directory

        Raises
        ------
        FileExistsError
             Raised if it is not a directory.
        PermissionError
             Raised if not possible to create the directory.
        """
        if self.data_dir is None:
            self.data_dir = self.default_data_dir

        if not self.data_dir.is_dir():
            try:
                self.data_dir.mkdir(parents=True)
            except FileExistsError as err:
                raise SpextraError(f"{self.data_dir} not a directory") from err
            except PermissionError as err:
                raise SpextraError(f"Cannot create {self.data_dir}") from err

        return self.data_dir

    def get_database_url(self) -> str:
        """
        Retrieve the database url and check whether is reachable.

        Returns
        -------
        url: str
            The url location of the database.
        """
        loc = self.database_url
        try:
            request = Request(loc)
            request.get_method = lambda: "HEAD"
            urlopen(request)
            is_url = True
        except URLError:
            is_url = False
        except ValueError:
            is_url = False

        if is_url is False:
            warnings.warn("Database might not be reachable")

        return loc

    def _set_param(self, name, value) -> None:
        dic = {name: value}
        self.meta.update(dic)
        self.__dict__.update(dic)
        with self.config_file.open("w", encoding="utf-8") as file:
            yaml.dump(self.meta, file, sort_keys=False)

    def __repr__(self) -> str:
        return yaml.dump(self.meta, indent=4, sort_keys=False,
                         default_flow_style=False)


def _download_file(remote_url, target, silent=False):
    """
    Accepts a URL, downloads the file to a given open file object.

    This is a modified version of astropy.utils.data.download_file that
    downloads to an open file object instead of a cache directory.
    """

    conf = Config()
    timeout = conf.remote_timeout
    download_block_size = 32768

    msg_file = None
    if silent:  # write the messages to dev/null
        msg_file = open(devnull, "w")

    try:
        # Pretend to be a web browser (IE 6.0). Some servers that we download
        # from forbid access from programs.
        headers = {"User-Agent": "Mozilla/5.0",
                   "Accept": ("text/html,application/xhtml+xml,"
                              "application/xml;q=0.9,*/*;q=0.8")}
        req = Request(remote_url, headers=headers)
        with closing(urlopen(req, timeout=timeout)) as remote:
            # get size of remote if available (for use in progress bar)
            info = remote.info()
            size = None
            if "Content-Length" in info:
                try:
                    size = int(info["Content-Length"])
                except ValueError:
                    pass

            dlmsg = f"Downloading {remote_url}"

            with ProgressBarOrSpinner(size, dlmsg, file=msg_file) as p:
                bytes_read = 0
                block = remote.read(download_block_size)
                while block:
                    target.write(block)
                    bytes_read += len(block)
                    p.update(bytes_read)
                    block = remote.read(download_block_size)

    # Append a more informative error message to HTTPErrors, URLErrors.
    except HTTPError as err:
        err.msg += f". requested URL: {remote_url!r}"
        raise err
    except URLError as err:
        append_msg = (hasattr(err, "reason") and hasattr(err.reason, "errno")
                      and err.reason.errno == 8)
        if append_msg:
            msg = f"{err.reason.strerror}. requested URL: {remote_url}"
            err.reason.strerror = msg
            err.reason.args = (err.reason.errno, msg)
        raise err

    # This isn't supposed to happen, but occasionally a socket.timeout gets
    # through.  It's supposed to be caught in `urrlib2` and raised in this
    # way, but for some reason in mysterious circumstances it doesn't. So
    # we'll just re-raise it here instead.
    except socket.timeout as err:
        # add the requested URL to the message (normally just "timed out")
        err.args = (f"requested URL {remote_url!r} timed out",)
        raise URLError(err)


def download_file(remote_url, local_name, silent=False):
    """
    Download a remote file to local path
    ----------
    remote_url : str
        The URL of the file to download
    local_name : str
        Absolute path filename of target file.
    Raises
    ------
    URLError
        Whenever there's a problem getting the remote file.
    """
    # ensure target directory exists
    path = Path(local_name)
    path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with path.open("wb") as target:
            _download_file(remote_url, target, silent=silent)
    except Exception as err:  # noqa
        # in case of error downloading, remove file.
        path.unlink(missing_ok=True)
        raise err


def download_svo_filter(filter_name):
    """
    Query the SVO service for the true transmittance for a given filter
    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the spanish VO filter service
        e.g: ``Paranal/HAWKI.Ks``

    Returns
    -------
    filt_curve : tuple with wave and trans values
    """
    from astropy.table import Table

    conf = Config()
    data_dir = conf.get_data_dir()
    data_dir_package = __data_dir__ / "svo"

    origin = ("http://svo2.cab.inta-csic.es/theory/fps3/"
              f"fps.php?ID={filter_name}")

    local_path_cache = Path(data_dir, "svo_filters", filter_name)
    local_path_package = Path(data_dir_package, filter_name)

    if local_path_package.exists():
        local_path = local_path_package
    else:
        local_path = local_path_cache
        if not local_path.exists():
            download_file(origin, local_path)

    # raises ValueError if table is malformed
    # this can be used to catch problmes
    tbl = Table.read(local_path, format="votable")

    wave = tbl["Wavelength"]
    trans = tbl["Transmission"]

    return wave, trans


def dict_generator(indict, pre=None):
    """
    Make a generator out of a dictionary.

    Parameters
    ----------
    indict: dict
    pre

    Yields
    ------
    list
    """
    pre = pre[:] if pre else []
    if isinstance(indict, dict):
        for key, value in indict.items():
            if isinstance(value, dict):
                yield from dict_generator(value, pre + [key])
            elif isinstance(value, (tuple, list)):
                for val in value:
                    yield from dict_generator(val, pre + [key])
            else:
                yield pre + [key]
    else:
        yield pre + [indict]


def get_filter_systems():
    """
    Return a set of the different filter system available.

    Returns
    -------

    """
    try:
        import tynt
    except ImportError as err:
        print(err, "this function requires tynt.")
    filters = tynt.FilterGenerator().available_filters()
    systems = {f.split("/")[0] for f in filters}
    return systems


def get_filter_names(system=None):
    """
    Return the filters available from tynt.

    if system= None returns all

    Returns
    -------

    """
    try:
        import tynt
    except ImportError as err:
        print(err, "this function requires tynt.")
    filter_list = tynt.FilterGenerator().available_filters()
    ord_list = [[f for f in filter_list if s in f]
                for s in get_filter_systems()]
    flat_list = [item for sublist in ord_list for item in sublist]

    if system is not None:
        flat_list = [f for f in filter_list if system in f]

    return flat_list


def _ensure_list(value):
    if isinstance(value, ARRAY_LIKE):
        return list(value)
    return [value]       


def _angstrom_qty(value):
    if isinstance(value, u.Quantity):
        return value.to(u.AA, equivalencies=u.spectral())
    return value * u.AA


def _angstrom_value(value):
    if isinstance(value, u.Quantity):
        return _angstrom_qty(value).value
    return value


def _abmag_qty(value):
    if not isinstance(value, u.Quantity):
        return value * u.ABmag
    return value
