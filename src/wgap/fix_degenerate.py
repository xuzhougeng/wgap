#!/usr/bin/env python

# MAKER only support A, C, G, T, N
# This script will convert all other characters(RYSWKMBDHV) to N

import sys

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: fix_degenerate.py <input fasta file> <output fasta file>")
        sys.exit(1)

    fasta = sys.argv[1]
    out_file = sys.argv[2]

    base = ['A', 'T', 'C', 'G', 'N']

    with open(out_file, 'w' ) as handler:
        for line in open(fasta, 'r') :
            if line.startswith('>'):
                handler.write(line)
                continue

            new_line = ''
            for char in line.strip():
                if char.upper() not in  base:
                    new_line += 'N'
                else:
                    new_line += char
            handler.write(new_line + '\n')