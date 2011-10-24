#!/usr/bin/env python

# based on Jeet Sukumaran's code
# https://gist.github.com/1015874

import os
import sys
import dendropy
import argparse

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_args():
    parser = argparse.ArgumentParser(description="""splits a nexus file by character
        partition""")
    parser.add_argument("nexus", help="nexus file to split", action=FullPaths)
    parser.add_argument("--output", help="folder to write split files to", action=FullPaths, type=is_dir)
    return parser.parse_args()

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass

def main():
    args = get_args()
    full_data = dendropy.DnaCharacterMatrix.get_from_path(args.nexus, 'nexus')
    items = full_data.character_subsets.items()
    (input_folder, input_name) = os.path.split(args.nexus)
    if args.output:
        output_folder = args.output
    else:
        output_folder = input_folder
    basename = os.path.splitext(input_name)[0]
    print "Writing to {0}".format(output_folder)
    mkdir(output_folder)
    for idx, (char_subset_name, char_subset) in enumerate(items):
        print "{0}/{1}: {2}".format(idx+1, len(items), char_subset_name)
        subset_data = full_data.export_character_subset(char_subset)
        output_name = "{0}.{1}.nex".format(basename, char_subset_name)
        output_full_path = os.path.join(output_folder, output_name)
        subset_data.write_to_path(output_full_path, 'nexus')

if __name__ == "__main__":
    main()

