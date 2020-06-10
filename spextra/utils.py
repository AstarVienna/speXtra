import os
import inspect
from astropy.config import ConfigItem, ConfigNamespace, get_cache_dir
from contextlib import closing
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from astropy.utils.console import ProgressBarOrSpinner


__all__ = ["is_url", "database_url", "download_file",  "Conf", "get_rootdir"]

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")


def is_url(url):  #
    """
    Checks that a given URL is reachable.
    Depending of the configuration of the server it might return True even if the page doesn't exist.

    Parameters
    -----------
    url: A URL

    Returns
    -------
    Boolean
    """
    try:
        request = Request(url)
        request.get_method = lambda: 'HEAD'
        urlopen(request)
        output = True
    except URLError:
        output = False
    except ValueError:
        output = False
    finally:
        return output


def database_url():
    """
    TODO: it should read it from a file
    Returns
    -------
    the database_location
    """
    loc = "https://homepage.univie.ac.at/miguel.verdugo/database/"
    try:
        assert is_url(loc)
    except AssertionError as error:
        print(error)
        print("Database address not reachable", loc)
    return loc


# ------ shamefully copied from SNCOSMO ------
class Conf(ConfigNamespace):
    """Configuration parameters for spextra."""
    data_dir = ConfigItem(
        None,
        "Directory where spextra will store and read downloaded data "
        "resources. If None, ASTROPY_CACHE_DIR/spextra is created and "
        "used. Example: data_dir = /home/user/data/sncosmo",
        cfgtype='string(default=None)')

    remote_timeout = ConfigItem(
        10.0, "Remote timeout in seconds.")


# Create an instance of the class we just defined.
# This needs to be done before the imports below because `conf` is used
# in some parts of the library.

conf = Conf()



def get_rootdir():
    # use the environment variable if set
    data_dir = os.environ.get('SPEXTRA_DATA_DIR')

    # otherwise, use config file value if set.

    if data_dir is None:
        data_dir = conf.data_dir

    # if still None, use astropy cache dir (and create if necessary!)
    if data_dir is None:
        data_dir = os.path.join(get_cache_dir(), "spextra")
        if not os.path.isdir(data_dir):
            if os.path.exists(data_dir):
                raise RuntimeError("{0} not a directory".format(data_dir))
            os.mkdir(data_dir)

    return data_dir


def _download_file(remote_url, target):
    """
    Accepts a URL, downloads the file to a given open file object.
    This is a modified version of astropy.utils.data.download_file that
    downloads to an open file object instead of a cache directory.
    """

    timeout = conf.remote_timeout
    download_block_size = 32768
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
            with ProgressBarOrSpinner(size, dlmsg) as p:
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


def download_file(remote_url, local_name):
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
    if not os.path.exists(dn):
        os.makedirs(dn)

    try:
        with open(local_name, 'wb') as target:
            _download_file(remote_url, target)
    except:  # noqa
        # in case of error downloading, remove file.
        if os.path.exists(local_name):
            os.remove(local_name)
            raise
