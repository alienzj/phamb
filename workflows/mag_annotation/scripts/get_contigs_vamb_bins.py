#!/usr/bin/python
import sys
import argparse
import vambtools as _vambtools
import run_RF_modules
import collections as _collections
import os
import numpy as _np


parser = argparse.ArgumentParser(
    description="""Command-line benchmark utility.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True)

parser.add_argument('fastafile', help='Path to concatenated assembly VAMB')
parser.add_argument('clusterspath', help='Path to clusters.tsv')
parser.add_argument('directoryout', help='Path to directory out', default="resultdir")
parser.add_argument('-m', dest='min_bin_size', metavar='', type=int, default=5000, help='Minimum size of bins - default [5000]')
parser.add_argument('-b', dest='batch_num', type=int, default=1000, help='number of contigs for each fna file' )


def write_concat_bins(directory, bins, fastadict, compressed=False, maxbins=250, minsize=5000, batch_num=1000):
    """Writes bins as FASTA files in a directory, one file per bin.
    Inputs:
        directory: Directory to create or put files in
        bins: {'name': {set of contignames}} dictionary (can be loaded from
        clusters.tsv using vamb.cluster.read_clusters)
        fastadict: {contigname: FastaEntry} dict as made by `loadfasta`
        compressed: Sequences in dict are compressed [False]
        maxbins: None or else raise an error if trying to make more bins than this [250]
        minsize: Minimum number of nucleotides in cluster to be output [0]
    Output: None
    """
    import os as _os
    import gzip as _gzip
    import vambtools as _vambtools

    import random

    # Safety measure so someone doesn't accidentally make 50000 tiny bins
    # If you do this on a compute cluster it can grind the entire cluster to
    # a halt and piss people off like you wouldn't believe.
    if maxbins is not None and len(bins) > maxbins:
        raise ValueError('{} bins exceed maxbins of {}'.format(len(bins), maxbins))

    # Check that the directory is not a non-directory file,
    # and that its parent directory indeed exists
    abspath = _os.path.abspath(directory)
    parentdir = _os.path.dirname(abspath)

    if parentdir != '' and not _os.path.isdir(parentdir):
        raise NotADirectoryError(parentdir)

    if _os.path.isfile(abspath):
        raise NotADirectoryError(abspath)

    if minsize < 0:
        raise ValueError("Minsize must be nonnegative")

    # Check that all contigs in all bins are in the fastadict
    allcontigs = set()

    for contigs in bins.values():
        allcontigs.update(set(contigs))

    allcontigs -= fastadict.keys()
    if allcontigs:
        nmissing = len(allcontigs)
        raise IndexError('{} contigs in bins missing from fastadict'.format(nmissing))

    # Make the directory if it does not exist - if it does, do nothing
    try:
        _os.mkdir(directory)
    except FileExistsError:
        pass

    bins_entries = []
    # Now actually print all the contigs to files
    for binname, contigs in bins.items():

        # Concatenate sequences of the bin
        concat_sequence = bytearray()
        for contig in contigs:
            entry = fastadict[contig]
            if compressed:
                uncompressed = bytearray(_gzip.decompress(entry.sequence))
                concat_sequence += uncompressed
            else:
                uncompressed = bytearray(entry.sequence)
                concat_sequence += uncompressed

        bin_entry = _vambtools.FastaEntry(binname, concat_sequence)
        # Skip bin if it's too small
        if len(bin_entry.sequence) < minsize:
            continue
        bins_entries.append(bin_entry)

    random.shuffle(bins_entries)
    print('Writing:',len(bins_entries) ,'bins to file')
    filename = _os.path.join(directory, 'vamb_bins.1.fna')
    i = 1
    j = 1
    file = open(filename,'w')
    for entry in bins_entries:
        if i % batch_num == 0:
            j += 1
            file.close()
            filename = _os.path.join(directory, 'vamb_bins.' + str(j) + '.fna')
            file = open(filename,'w')
        i += 1
        print(entry.format(),file=file)


if __name__=='__main__':

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    try:
        os.mkdir(args.directoryout)
    except FileExistsError:
        pass

    with open(args.clusterspath) as file:
        clusters = _vambtools.read_clusters(file)

    with _vambtools.Reader(args.fastafile, 'rb') as infile:
        fastadict = _vambtools.loadfasta(infile,compress=False)

    reference = run_RF_modules.Reference.from_clusters(
        clusters = clusters,
        fastadict=fastadict,
        minimum_contig_len=2000)

    bins_all = {binname:clusters[binname] for binname in reference.genomes.keys()}
    print(len(bins_all))

    write_concat_bins(
        args.directoryout,
        bins_all, fastadict,
        compressed=False,
        maxbins=len(bins_all),
        minsize=args.min_bin_size,
        batch_num=args.batch_num)
