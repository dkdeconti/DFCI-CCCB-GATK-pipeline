#! /usr/bin/python

import re
import sys


def fix_line(line):
    arow = line.strip('\n').split('\t')
    sample_name = arow[3][3:]
    sample_name = sample_name[:re.search("_S", sample_name).start()]
    arow[3] = "LB:" + sample_name
    arow[4] = "SM:" + sample_name
    return '\t'.join(arow) + '\n'


def main(sa):
    header_filename = sa[0]
    with open(header_filename, 'rU') as handle:
        for line in handle:
            if line[0:3] == "@RG":
                sys.stdout.write(fix_line(line))
            else:
                sys.stdout.write(line)


if __name__ == "__main__":
    main(sys.argv[1:])

