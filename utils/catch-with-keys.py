# Copyright (C) 2017 Zhixian MA <zxma_sjtu@qq.com>

"""
Catch data with provided keys

Reference
=========
http://www.cnblogs.com/jiu0821/p/6275685.html
"""

import os
import re
from pandas import DataFrame

def main():
    keys = [r"n0", r"beta", r"rc_kpc", r"rho0", r"rs"]

    # load file
    filename = "./mass_20170909.log"
    print("Loading the file %s" % (filename))
    fp = open(filename)
    lines = fp.readlines()

    print("Catching lines with keys ", keys)
    # catch data
    datadict = {"n0":[],"beta":[],"rc_kpc":[],"rho0":[],"rs":[]}
    for l in lines:
        for k in keys:
            d = re.findall(k, l)
            if len(d):
                if len(l.split('\t')) == 2:
                    datadict[k].append(float(l.split('\t')[1]))

    fp.close()

    # Do not save the last pair
    del(datadict["beta"][-1])
    del(datadict["n0"][-1])
    del(datadict["rc_kpc"][-1])

    # save result
    savename = filename[0:-4] + "_catch.csv"
    print("Saving result to %s" % (savename))
    dataframe = DataFrame(datadict)
    dataframe.to_csv(savename, sep=" ")

if __name__ == "__main__":
    main()
