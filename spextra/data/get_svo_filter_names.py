"""
script to construct the list of filters available at SVO
"""
import os
import inspect

from spextra.utils import get_rootdir,  download_file

from astropy.table import Table

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# Default filters


def get_names():
    facilities_list = []
    with open("speXtra/spextra/data/facilities_list.txt") as facilities:
        for line in facilities:
            facilities_list.append(line[:-1])

    print(facilities_list)
    filter_IDs = []
    filename = os.path.join(get_rootdir(), "temp_table" )
    for f in facilities_list:
        search_form = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?Facility={}".format(f)
        table_file = download_file(search_form, filename)
        try:
            table = Table.read(filename) 
            filter_IDs.extend(list(table["filterID"]))
        except ValueError:
            print(f)
            pass

    return filter_IDs


def main():
    filter_names = get_names()
    with open("filter_names.txt", "w") as f:
        for item in filter_names:
            f.write('%s\n' % item)


if __name__ == "__main__":
    main()
