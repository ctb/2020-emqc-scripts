#! /usr/bin/env python
import sourmash
import sys
import argparse
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs')   # this is an assembly
    p.add_argument('read_sig')  #  this contains sourmash sig with abunds
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    siglist = sourmash.load_file_as_signatures(args.read_sig)
    siglist = list(siglist)
    assert len(siglist) == 1
    sig = siglist[0]

    contigs_mh = sig.minhash.copy_and_clear()
    for record in screed.open(args.contigs):
        contigs_mh.add_sequence(record.sequence, force=True)

    # intersect the genome assembly with the read abundances
    # so now we get the abundances of only the k-mers that are in the
    # assembly.
    abunds = {}
    for hashval in contigs_mh.hashes:
        abunds[hashval] = sig.minhash.hashes.get(hashval, 0)

    output_mh = sig.minhash.copy_and_clear()
    output_mh.set_abundances(abunds)

    out_sig = sourmash.SourmashSignature(output_mh)
    with open(args.output, 'wt') as fp:
        print(f"Saving output to '{args.output}'")
        sourmash.save_signatures([out_sig], fp)


if __name__ == '__main__':
    main()
