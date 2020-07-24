"""
script to construct the list of filters available at SVO
"""
import os
import inspect

from .utils import get_rootdir,  download_file

from astropy.table import Table

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# Default filters


def get_names():
    with open(os.path.join(__data_dir__,  "facilities_list.txt")) as facilities:
        facilities_list = facilities.read()

    filter_IDs = []
    for f in facilities_list:
        search_form = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?Facility={}".format(f)
        table_file = download_file(search_form, get_rootdir())
        table = Table.read(table_file)
        filter_IDs.extend(list(table["FilterID"]))

    return filter_IDs


def main():
    filter_names = get_names()
    with open(os.path.join(__data_dir__,  "filter_names.txt"), "w") as f:
        for item in filter_names:
            f.write('%s\n' % item)


if __name__ == "__main__":
    main()
