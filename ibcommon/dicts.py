def merge_dicts(dicts):
    """Given list of dicts, merge them into a new dict as a shallow copy.
    :param dicts: list of dictionaries
    :return: merged dictionary
    """
    out_dict = {}
    for current_dict in dicts:
        if not isinstance(current_dict,dict):
            continue
        out_dict = out_dict.copy()
        try:
            out_dict.update(current_dict)
        except:
            print 'error append dict'
    return out_dict


def merge_two_dicts(dict_1, dict_2):
    """Given two dicts, merge them into a new dict as a shallow copy.
    :param dict_1: dictionary 1
    :param dict_2: dictionary 2
    :return: merged dictionary
    """
    if not dict_1:
        return dict_2
    if not dict_2:
        return dict_1

    out_dict = dict_1.copy()
    out_dict.update(dict_2)
    return out_dict


def extract_keys(inp_dict, keys):
    """
    Return dictionary witxh selected keys
    :param inp_dict: dictionary
    :param keys: list or string
    :return: dictionary
    """
    if not isinstance(keys, list):
        keys = [keys]
    return {your_key: inp_dict[your_key] for your_key in keys}


def get_column(inp_dict, key):
    """
    return one 'column' from array of dictionaries
    :param inp_dict: dictionary 
    :param key: column key
    :return: array
    """
    return [row[key] for row in inp_dict]
