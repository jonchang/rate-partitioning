#!/usr/bin/env python

import os
import csv
import operator
import argparse
import collections
import itertools

import dendropy
import numpy
from scipy.cluster.vq import kmeans2

import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """Get arguments / options"""
    parser = argparse.ArgumentParser(description="""creates nexus file
        partitions based on rates""")
    parser.add_argument('sites', help="some file with sequences",
        action=FullPaths)
    parser.add_argument('rates', help="""rates in a tab-separated or
        comma-separated file, where the first field is the site number,
        and the second field is the rate for that site""", action=FullPaths)
    parser.add_argument('--site-format', help="format for the sites file",
        default="nexus", choices=["nexus"])
    parser.add_argument('--clusters', help="""number of clusters to partition
        the data into""", required=True, type=int)
    parser.add_argument('--output', help="""name of nexus file to output the
        clustered rates""", action=FullPaths, default="output.nex")
    parser.add_argument("--rearrange-sites", help="""rearrange the sites in the
        output nexus file so that clusters are contiguous partitions""",
        action="store_true")

    return parser.parse_args()

def guess_csv(filehandle):
    """Sniffs a csv file's format and guesses if it has a header or not"""
    sniffer = csv.Sniffer()
    sample = filehandle.read(1024)
    filehandle.seek(0)
    dialect = sniffer.sniff(sample, delimiters="\t,")
    header = sniffer.has_header(sample)
    return dialect, header

def group_ranges(L):
    """groups a list L into a list of ranges where each range has integers
    that differ from the previous integer in that range by 1"""
    for w, z in itertools.groupby(L, lambda x, y=itertools.count(): next(y) - x):
        grouped = list(z)
        yield [str(x) for x in [grouped[0], grouped[-1]][:len(grouped)]]

def range_to_string(L):
    """converts a list of ranges from group_ranges into a string suitable for
    consumption by the nexus format."""
    groups = list(group_ranges(L))
    return ' '.join(['-'.join(x) for x in groups])

def create_sets_block(sets):
    """Creates a sets block. `sets` is a dictionary where keys are character
    set names and values are a list of character indicies for that set."""
    lines = []
    lines.append("begin sets;")
    for k, v in sets.iteritems():
        lines.append("charset cluster_{0} = {1};".format(k, range_to_string(v)))
    lines.append("end;")
    return "\n".join(lines)

def cluster(values, num_clusters):
    """Turn `values` into `num_clusters` clusters"""
    value_range = max(values) - min(values)
    posts = [min(values)]
    for idx in range(1, num_clusters):
        posts.append(float(idx) / num_clusters * value_range)
    posts.append(max(values))

    obs = []
    for idx, value in enumerate(posts, start=1):
        if idx < len(posts):
            obs.append(value + posts[idx] / 2)

    centroids, labels = kmeans2(numpy.array(values), numpy.array(obs), minit="matrix")
    clusters = collections.defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(idx + 1) # sites start counting at 1

    return clusters

def get_rates(filename):
    """Gets site rates from a rate file `filename` and returns it as a dict"""
    with open(filename, "rU") as rfile:
        dialect, header = guess_csv(rfile)
        csvreader = csv.reader(rfile, dialect=dialect)
        if header: # skip header row
            csvreader.next()
        rates = {}
        for row in csvreader:
            rates[row[0]] = float(row[1])
        return rates

def main():
    args = get_args()

    all_rates = get_rates(args.rates)

    # Sort rates correctly
    sorted_rates = collections.OrderedDict(sorted(all_rates.iteritems(),
                                                  key=lambda x: int(x[0])))

    # Pull out values in the correct order
    clusters = cluster(sorted_rates.values(), args.clusters)

    # Write output data
    data = dendropy.DataSet.get_from_path(args.sites, 'nexus')
    if args.rearrange_sites:
        matrix = data.char_matrices[0]
        # python lists start at 0
        cluster_values = []
        for x in clusters.itervalues():
            for y in x:
                cluster_values.append(y - 1)
        getters = [operator.itemgetter(x) for x in cluster_values]
        print '\n'.join([str(max(x)) for x in clusters.itervalues()])
        for taxon, characters in matrix.iteritems():
            new_order = [x(characters) for x in getters]
            pdb.set_trace()
    data.write_to_path(args.output, "nexus", supplemental_blocks=[create_sets_block(clusters)])

if __name__ == "__main__":
    main()
