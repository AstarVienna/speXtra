import tynt
from spextra.utils import download_file

def get_filter_systems():
    """
    Return a set of the different filter system available

    Returns
    -------

    """
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
    filter_list = tynt.FilterGenerator().available_filters()
    ord_list = [[f for f in filter_list if s in f] for s in get_filter_systems()]
    flat_list = [item for sublist in ord_list for item in sublist]

    if system is not None:
        flat_list = [f for f in filter_list if system in f]

    return flat_list


def download_filter(filter_name, path):
    path = download_file('http://svo2.cab.inta-csic.es/'
                          'theory/fps3/fps.php?ID={}'.format(filter_name), path)
    


if __name__ == '__main__':
    system_list = get_filter_systems()
    for s in system_list:
        filter_list = get_filter_names(system=s)
        for f in filter_list:
            download_filter(filter_name=f, path=f+".vo")



    
    


    


