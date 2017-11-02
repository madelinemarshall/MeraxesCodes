#!/usr/bin/env python3

"""Get the on switches for a merger tree flag integer.

Usage: parse_tree_flags.py <flag>
"""

from docopt import docopt
import sys

__author__ = "Simon Mutch"
__date__ = "21 Jan 2014"


class TreeFlags():

    def __init__(self):
        self.flags = {}
        self.values = {}
        with open("tree_flags.h", "r") as fd:
            for line in fd:
                line = line.split()
                if len(line) == 0:
                    continue
                if not line[0].startswith('//'):
                    if line[1].startswith('TTT'):
                        self.values[line[1]] = int(line[2])
                    else:
                        self.flags[line[1]] = line[2]

        for k, v in self.flags.items():
            if isinstance(v, str):
                self.flags[k] = self.values[v]

    def parse(self, num):
        print("Flag {:d} matches:".format(num))
        for s, v in self.flags.items():
            if (num & v) == v:
                print(s)
   
    def flaglist(self, num):
        flaglist=[]
        for s, v in self.flags.items():
            if (num & v) ==v:
                flaglist.append(str(s))
        return(flaglist)

if __name__ == '__main__':
    #args = docopt(__doc__)
    flags = TreeFlags()
    flags.parse(int(sys.argv[1]))
