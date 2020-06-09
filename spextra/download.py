import shutil
from posixpath import join as urljoin
import urllib
import os
import inspect
from astropy.config import ConfigItem, ConfigNamespace, get_cache_dir

#from astropy.utils.data import download_file, clear_download_cache, 
from astropy.utils.data import clear_download_cache, cache_contents
import yaml
#clear_download_cache()


# %%%%% shamefully copied from SNCOSMO
class _Conf(ConfigNamespace):
    """Configuration parameters for sncosmo."""
    data_dir = ConfigItem(
        None,
        "Directory where sncosmo will store and read downloaded data "
        "resources. If None, ASTROPY_CACHE_DIR/sncosmo is created and "
        "used. Example: data_dir = /home/user/data/sncosmo",
        cfgtype='string(default=None)')
    
    remote_timeout = ConfigItem(
        10.0, "Remote timeout in seconds.")


# Create an instance of the class we just defined.
# This needs to be done before the imports below because `conf` is used
# in some parts of the library.
conf = _Conf()


def get_rootdir():
    # use the environment variable if set
    data_dir = os.environ.get('SNCOSMO_DATA_DIR')

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

    from contextlib import closing
    from urllib.request import urlopen, Request
    from urllib.error import URLError, HTTPError
    from astropy.utils.console import ProgressBarOrSpinner
    from . import conf

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
    Download a remote file to local path, unzipping if the URL ends in '.gz'.
    Parameters
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

    if remote_url.endswith(".gz"):
        import io
        import gzip

        buf = io.BytesIO()
        _download_file(remote_url, buf)
        buf.seek(0)
        f = gzip.GzipFile(fileobj=buf, mode='rb')

        with open(local_name, 'wb') as target:
            target.write(f.read())
        f.close()

    else:
        try:
            with open(local_name, 'wb') as target:
                _download_file(remote_url, target)
        except:  # noqa
            # in case of error downloading, remove file.
            if os.path.exists(local_name):
                os.remove(local_name)
            raise


def download_dir(remote_url, dirname):
    """
    Download a remote tar file to a local directory.
    Parameters
    ----------
    remote_url : str
        The URL of the file to download
    dirname : str
        Directory in which to place contents of tarfile. Created if it
        doesn't exist.
    Raises
    ------
    URLError (from urllib2 on PY2, urllib.request on PY3)
        Whenever there's a problem getting the remote file.
    """

    import io
    import tarfile

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    mode = 'r:gz' if remote_url.endswith(".gz") else None

    # download file to buffer
    buf = io.BytesIO()
    _download_file(remote_url, buf)
    buf.seek(0)

    # create a tarfile with the buffer and extract
    tf = tarfile.open(fileobj=buf, mode=mode)
    tf.extractall(path=dirname)
    tf.close()
    buf.close()  # buf not closed when tf is closed.


class DataMirror(object):
    """
    TODO: merge in SpecDatabase
    
    Lazy fetcher for remote data.
    When asked for local absolute path to a file or directory, DataMirror
    checks if the file or directory exists locally and, if so, returns it.
    If it doesn't exist, it first determines where to get it from.
    It first downloads the file ``{remote_root}/redirects.json`` and checks
    it for a redirect from ``{relative_path}`` to a full URL. If no redirect
    exists, it uses ``{remote_root}/{relative_path}`` as the URL.
    It downloads then downloads the URL to ``{rootdir}/{relative_path}``.
    For directories, ``.tar.gz`` is appended to the
    ``{relative_path}`` before the above is done and then the
    directory is unpacked locally.
    Parameters
    ----------
    rootdir : str or callable
        The local root directory, or a callable that returns the local root
        directory given no parameters. (The result of the call is cached.)
        Using a callable allows one to customize the discovery of the root
        directory (e.g., from a config file), and to defer that discovery
        until it is needed.
    remote_root : str
        Root URL of the remote server.
    """

    def __init__(self, rootdir, remote_root):
        if not remote_root.endswith('/'):
            remote_root = remote_root + '/'

        self._checked_rootdir = None
        self.rootdir = rootdir
        self.remote_root = remote_root

    

    def rootdir(self):
        """Return the path to the local data directory, ensuring that it
        exists"""

        if self._checked_rootdir is None:

            # If the supplied value is a string, use it. Otherwise
            # assume it is a callable that returns a string)
            rootdir = (self._rootdir
                       if isinstance(self._rootdir, str)
                       else self._rootdir())

            # Check existance
            if not os.path.isdir(rootdir):
                raise Exception("data directory {!r} not an existing "
                                "directory".format(rootdir))

            # Cache value for future calls
            self._checked_rootdir = rootdir

        return self._checked_rootdir



    def abspath(self, relpath, isdir=False):
        """Return absolute path to file or directory, ensuring that it exists.
        If ``isdir``, look for ``{relpath}.tar.gz`` on the remote server and
        unpackage it.
        Otherwise, just look for ``{relpath}``. If redirect points to a gz, it
        will be uncompressed."""

        abspath = os.path.join(self.rootdir(), relpath)

        if not os.path.exists(abspath):
            if isdir:
                url = urljoin(self.remote_root, relpath)
                
                # Download and unpack a directory.
                download_dir(url, os.path.dirname(abspath))

                # ensure that tarfile unpacked into the expected directory
                if not os.path.exists(abspath):
                    raise RuntimeError("Tarfile not unpacked into expected "
                                       "subdirectory. Please file an issue.")
            else:
                url = urljoin(self.remote_root, relpath)
                download_file(url, abspath)

        return abspath

#

rootdir = get_rootdir()
DATADIR = DataMirror(get_rootdir, "https://homepage.univie.ac.at/miguel.verdugo/database/")


def get_yaml_contents(relpath):
    """
    read a yaml file from a relative url

    Parameters
    ----------
    path

    Returns
    -------
    dict with the contents of the yaml file
    """


    filename = DATADIR.abspath(relpath)
    with open(filename) as f:
        data = yaml.safe_load(f)

    return data


def get_library(library_name):
    
    relpath = urljoin("libraries", library_name, library_name + ".yml")

    data = get_yaml_contents(relpath)
    
    return data

def get_template(template):
    """
    

    Parameters
    ----------
    relpath : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    library_name, template_name = template.split("/")
    lib_data = database.get_library(library_name)
    template_meta = {"resolution": lib_data["resolution"],
                     "wave_unit": lib_data["wave_unit"],
                     "flux_unit": lib_data["flux_unit"],
                     "wave_column_name": lib_data["wave_column_name"],
                     "flux_column_name": lib_data["flux_column_name"],
                     "data_type": lib_data["data_type"],
                     "file_extension": lib_data["file_extension"]}
 
    filename = template_name + template_meta["file_extension"]
    relpath =  urljoin("libraries", library_name, filename)

    try:
        assert template_name in lib_data["templates"]
    except AssertionError as error:
        print(error)
        print(template_name, "not found")
    else:
        filename = DATADIR.abspath(relpath)

    return newfile, template_meta



print(rootdir)



print("*"*15)
print(DATADIR.rootdir())
print("*"*15)
print(DATADIR.remote_root)

print("%"*20)

#abspath = DATADIR.abspath("libraries/kc96/index.yml")
#print(abspath)

data = get_yaml_contents("libraries/kc96/index.yml")
print(data)
#remote_file = "https://homepage.univie.ac.at/miguel.verdugo/database/libraries/kc96/index.yml"
#file = download_file(remote_file, rootdir)

#print(file)

print(cache_contents("spextra"), "*"*15)


