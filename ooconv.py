import openpyxl as xl
import datetime
import string
import re

from tkinter import Tk
from tkinter.filedialog import askopenfilenames

header = 'headers.oo'

def readrow(ws, cols, i):
    row = []
    for letter in cols:
        val = ws[letter + str(i)].value
        row.append(val)
    return(row)

def xl2list(xlfile):
    # starting values
    start = 3
    oocols = string.ascii_uppercase
    oononelist = [None] * 26
    oorow = ''
    expcols = ['AA','AB','AC','AD','AE','AF','AG','AH','AI','AJ']
    expnonelist = [None] * 10
    exprow = ''
    fits = [None] * 6
    prefix = ['OOCYTE']
    # end starting values
    wb = xl.load_workbook(xlfile)
    ws = wb['DataEntry']
    i = start
    oolist = []
    while oorow != oononelist:
        oorow = readrow(ws, oocols, i)
        oolist.append(oorow)
        i += 1
    oolist.pop(-1)
    i = start
    explist = []
    while exprow != expnonelist:
        exprow = readrow(ws, expcols, i)
        explist.append(exprow)
        i += 1
    explist.pop(-1)
    dlist = []
    if len(explist) == 1:
        for row in oolist:
            newrow = prefix + row + fits + explist[0]
            dlist.append(newrow)
    elif len(explist) == len(oolist):
        for i in range(0, len(oolist)):
            newrow = prefix + row[i] + fits + explist[i]
            dlist.append(newrow)
    else:
        explist = explist[0]
        for row in oolist:
            newrow = prefix + row + fits + explist
            dlist.append(newrow)
    return(dlist)
def writeoofile(dlist, header, xlfile):
    with open(header, 'r') as hfile:
        header = hfile.readlines()
    hfile.close()
    oofile = dlist[0][1]
    p = re.compile(".+-.+-.+-.+-")
    res = p.findall(oofile)
    oofile = './dir/' + res[0][:-1] + '.oo'
    #print(oofile)
    with open(oofile, 'w') as opfile:
        for row in header:
            opfile.write(row)
        newdlist = []
        for row in dlist:
            newrow = []
            for value in row:
                if type(value) is int:
                    newval = str(value)
                elif type(value) is float:
                    newval = str(value)
                elif value is None:
                    newval = ""
                elif isinstance(value, datetime.datetime):
                    newval = value.strftime("%Y-%m-%d")
                else:
                    newval = value
                newrow.append(newval)
            newdlist.append(newrow)
            opfile.write(",".join(newrow) + "\n")
    opfile.close()
    print(oofile + ' written to disk')
    return(None)
def ooconv(xlfile, header):
    dlist = xl2list(xlfile)
    writeoofile(dlist, header, xlfile)
    return(None)
if __name__ == "__main__":
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    xlfiles = askopenfilenames(initialdir = ".",title = "Select a file to convert to an .oo", filetypes=[('Excel', ('*.xls', '*.xlsx'))])
    XLlist = list(xlfiles)
    for xlfile in XLlist:
        print(xlfile)
        ooconv(xlfile, header)
    

