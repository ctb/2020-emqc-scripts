#! /usr/bin/env python
import sys
import argparse
import sourmash
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sigfile')
    p.add_argument('contigs')
    p.add_argument('outfile')
    args = p.parse_args()

    sig = sourmash.load_one_signature(args.sigfile)
    print(f"loaded sig '{sig}' from {args.sigfile}")

    query_mh = sig.minhash

    print(f"saving matching contigs to {args.outfile}")
    outfp = open(args.outfile, 'wt')

    print(f"loading contigs from {args.contigs}")
    m = 0
    m_bp = 0
    for n, record in enumerate(screed.open(args.contigs)):
        if n and n % 1000 == 0:
            print(f'... {n + 1} contigs {m} saved ({m_bp/1000:.0f} kb)', end='\r', file=sys.stderr)

        new_mh = query_mh.copy_and_clear()
        new_mh.add_sequence(record.sequence.upper())
        if new_mh.count_common(query_mh) >= 5:
            m += 1
            m_bp += len(record.sequence)

            print(f">{record.name}\n{record.sequence}", file=outfp)

    print(f'done! {n + 1} contigs {m} saved ({m_bp/1000:.0f} kb)', file=sys.stderr)            

    return 0


if __name__ == '__main__':
    sys.exit(main())
