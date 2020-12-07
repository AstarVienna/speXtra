import os
import inspect
import warnings
from contextlib import closing
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import socket

from astropy.utils.console import ProgressBarOrSpinner
from astropy.config import get_cache_dir
import yaml


__all__ = ["Config", "download_file", "dict_generator", "download_svo_filter"]

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")
__config_file__ = os.path.join(__pkg_dir__, "data", "config.yml")


class Config:
    """
    set up configuration
    """

    def __init__(self, data_dir=None, database_url=None, remote_timeout=None):

        self.config_file = __config_file__

        with open(self.config_file) as f:
            self.meta = yaml.safe_load(f)
            self.data_dir = self.meta["data_dir"]
            self.database_url = self.meta["database_url"]
            self.remote_timeout = self.meta["remote_timeout"]
            self.default_data_dir = os.path.join(get_cache_dir(), self.meta["default_data_dir"])

        if data_dir is not None:
            self.set_param("data_dir", data_dir)

        if database_url is not None:
            self.set_param("database_url", database_url)

        if remote_timeout is not None:
            self.set_param("remote_timeout", remote_timeout)

    def get_data_dir(self):
        # use the environment variable if set
        data_dir = os.environ.get('SPEXTRA_DATA_DIR')

        # otherwise, use config file value if set.

        if data_dir is None:
            data_dir = self.data_dir

        # if still None, use astropy cache dir (and create if necessary!)
        if data_dir is None:
            data_dir = self.default_data_dir
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
            request = Request(loc)
            request.get_method = lambda: 'HEAD'
            urlopen(request)
            is_url = True
        except URLError:
            is_url = False
        except ValueError:
            is_url = False

        if is_url is False:
            warnings.warn("Database might not be reachable")

        return loc

    def set_param(self, name, value):
        d = {name: value}
        self.meta.update(d)
        self.__dict__.update(d)
        with open(self.config_file, "w") as f:
            yaml.dump(self.meta, f, sort_keys=False)

    def __repr__(self):
        return yaml.dump(self.meta, indent=4, sort_keys=False, default_flow_style=False)


def dict_generator(indict, pre=None):
    """
    Make a generator out of a dictionary
    Parameters
    ----------
    indict
    pre

    Returns
    -------

    """

    pre = pre[:] if pre else []
    if isinstance(indict, dict):
        for key, value in indict.items():
            if isinstance(value, dict):
                for d in dict_generator(value, pre + [key]):
                    yield d
            elif isinstance(value, list) or isinstance(value, tuple):
                for v in value:
                    for d in dict_generator(v, pre + [key]):
                        yield d
            else:
                yield pre + [key]
    else:
        yield pre + [indict]



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
    if silent is True:  # write the messages to dev/null
        msg_file = open(os.devnull, "w")

    try:
        # Pretend to be a web browser (IE 6.0). Some servers that we download
        # from forbid access from programs.
        headers = {'User-Agent': 'Mozilla/5.0',
                   'Accept': ('text/html,application/xhtml+xml,'
                              'application/xml;q=0.9,*/*;q=0.8')}
        req = Request(remote_url, headers=headers)
        with closing(urlopen(req, timeout=timeout)) as remote:

            # get size of remote if available (for use in progress bar)
            info = remote.info()
            size = None
            if 'Content-Length' in info:
                try:
                    size = int(info['Content-Length'])
                except ValueError:
                    pass

            dlmsg = "Downloading {0}".format(remote_url)

            with ProgressBarOrSpinner(size, dlmsg, file=msg_file) as p:
                bytes_read = 0
                block = remote.read(download_block_size)
                while block:
                    target.write(block)
                    bytes_read += len(block)
                    p.update(bytes_read)
                    block = remote.read(download_block_size)

    # Append a more informative error message to HTTPErrors, URLErrors.
    except HTTPError as e:
        e.msg = "{}. requested URL: {!r}".format(e.msg, remote_url)
        raise
    except URLError as e:
        append_msg = (hasattr(e, 'reason') and hasattr(e.reason, 'errno') and
                      e.reason.errno == 8)
        if append_msg:
            msg = "{0}. requested URL: {1}".format(e.reason.strerror,
                                                   remote_url)
            e.reason.strerror = msg
            e.reason.args = (e.reason.errno, msg)
        raise e

    # This isn't supposed to happen, but occasionally a socket.timeout gets
    # through.  It's supposed to be caught in `urrlib2` and raised in this
    # way, but for some reason in mysterious circumstances it doesn't. So
    # we'll just re-raise it here instead.
    except socket.timeout as e:
        # add the requested URL to the message (normally just 'timed out')
        e.args = ('requested URL {!r} timed out'.format(remote_url),)
        raise URLError(e)


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
    dn = os.path.dirname(local_name)
    if os.path.exists(dn) is False:
        os.makedirs(dn)

    try:
        with open(local_name, 'wb') as target:
            _download_file(remote_url, target, silent=silent)
    except:  # noqa
        # in case of error downloading, remove file.
        if os.path.exists(local_name):
            os.remove(local_name)
            raise


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

    origin = 'http://svo2.cab.inta-csic.es/'\
             'theory/fps3/fps.php?ID={}'.format(filter_name)

    local_path = os.path.join(conf.get_data_dir(), "svo_filters", filter_name)

    if os.path.exists(local_path) is False:
        download_file(origin, local_path)

    tbl = Table.read(local_path, format='votable')  # raises ValueError if table is malformed
                                                    # this can be used to catch problmes
    wave = tbl['Wavelength']
    trans = tbl['Transmission']

    return wave, trans


#    This is based on scopesim.effects.ter_curves_utils.py


def get_filter_systems():
    """
    Return a set of the different filter system available

    Returns
    -------

    """
    try:
        import tynt
    except ImportError as e:
        print(e, "this function requires tynt. ")
    filters = tynt.FilterGenerator().available_filters()
    systems = {f.split("/")[0] for f in filters}
    return systems


def get_filter_names(system=None):
    """
    This function just returns the filters available from tynt
    if system= None returns all

    Returns
    -------

    """
    try:
        import tynt
    except ImportError as e:
        print(e, "this function requires tynt. ")
    filter_list = tynt.FilterGenerator().available_filters()
    ord_list = [[f for f in filter_list if s in f] for s in get_filter_systems()]
    flat_list = [item for sublist in ord_list for item in sublist]

    if system is not None:
        flat_list = [f for f in filter_list if system in f]

    return flat_list