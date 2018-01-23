
import ConfigParser
import ast
import grp
import json
import os
import numpy as np

import xlrd


def parseExcel(filename, clmnnames=-1, datastart=0, sheetname='Sheet1', *argv):
    """
        Parse excel file into ....

        filename = excel file name
        clmnnames = row number where column names are stored (-1 if generic clmn1,clmn2...)
        datastart = row number where data starts
        sheetname = name of the sheet

    """
    try:
        book = xlrd.open_workbook(filename)
        sheet = book.sheet_by_name(sheetname)
    except:
        print "Problem opening " + filename
        return False
    types = sheet.row_types(datastart)
    if clmnnames == -1:
        nmcols = ['clmn' + str(k) for k in range(len(types))]
    else:
        nmcols = corrUniqCols([str(q) for q in sheet.row_values(clmnnames)])
    data = []
    for i in range(datastart, sheet.nrows):
        data.append(corrUnicodeErr(sheet.row_values(i)))
    return nmcols, types, data


def corrUniqCols(arr):
    newarr = []
    k = 0
    for cl in arr:
        k += 1
        newcl = cl.replace(" ", "_").lower()
        if newcl not in newarr and newcl != '':
            newarr.append(newcl)
        else:
            newarr.append('clmn_' + str(k) + "_" + newcl)

    return newarr


def corrUnicodeErr(row):
    newrow = []
    for q in row:
        try:
            str(q)
            newrow.append(q)
        except:
            newrow.append('-')
    return newrow


def pythonToPHP(result, tmp_file):
    to_php = json.dumps(result)
    file_down = open(tmp_file, 'w')
    file_down.write(to_php)
    file_down.close()
    gid = grp.getgrnam("brewers").gr_gid
    os.chown(tmp_file, -1, gid)
    os.chmod(tmp_file, 0777)


def PHPTopython(infile):
    return json.loads(open(infile).read(1000000))


def parse_ini_arguments(ini_file, section):
    config = ConfigParser.ConfigParser()
    config.read(ini_file)
    results = {arg: tryeval(config.get(section, arg)) for arg in config.options(section)}
    return results


def tryeval(val):
    """
    Convert to the proper type, e.g. '87' -> 87
    :param val: string
    :return: mixed
    """
    try:
        val = ast.literal_eval(val)
    except (SyntaxError, ValueError) as ex:
        pass
    return val


## correct if forward slash is missing at the end of 'path'
def corrpath(path):
    return os.path.join(path, '')


## parse True to 'y' and False to 'n'
def booltostring(inpbool):
    if inpbool:
        return 'y'
    return 'n'


## parse single vector to value i.e. [5] -> 5
def singlevec(value):
    if isinstance(value, np.ndarray) and len(value) == 1:
        return value[0]
    return value
