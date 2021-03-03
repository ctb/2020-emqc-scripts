#! /usr/bin/env python
"Extract contigs by their median k-mer abundance in reads sig."
import sys
import argparse
import sourmash
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sigfile')
    p.add_argument('contigs')
    p.add_argument('outfile')
    p.add_argument('--min', default=0, type=int)
    p.add_argument('--max', default=0, type=int)
    args = p.parse_args()

    sig = sourmash.load_one_signature(args.sigfile)
    print(f"loaded sig '{sig}' from {args.sigfile}")

    query_mh = sig.minhash
    assert query_mh.track_abundance, "reads sig must have abundance"

    print(f"saving contigs to {args.outfile}")
    outfp = open(args.outfile, 'wt')

    print(f"loading contigs from {args.contigs}")
    m = 0
    m_bp = 0
    for n, record in enumerate(screed.open(args.contigs)):
        if n and n % 1000 == 0:
            print(f'... {n + 1} contigs {m} saved ({m_bp/1000:.0f} kb)', end='\r', file=sys.stderr)

        new_mh = query_mh.copy_and_clear()
        new_mh.add_sequence(record.sequence.upper(), force=True)

        # for the k-mers in the contigs, get the abundances from the reads
        abunds = [ query_mh.hashes.get(x, 0) for x in new_mh.hashes ]
        abunds = list(sorted(abunds))

        if not abunds:
            continue

        median_abund = abunds[round(len(abunds) / 2)]
        if median_abund >= args.min and \
           (not args.max or median_abund <= args.max):
            m += 1
            m_bp += len(record.sequence)

            print(f">{record.name}\n{record.sequence}", file=outfp)

    print(f'done! {n + 1} contigs, {m} saved ({m_bp/1000:.0f} kb)', file=sys.stderr)            

    return 0


if __name__ == '__main__':
    sys.exit(main())
