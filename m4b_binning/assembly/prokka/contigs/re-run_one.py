import argparse

from prokka_set_of_files import run_prokka


if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    parser.add_argument('fasta', help='fasta to re-run')
    parser.add_argument('threads', help='number of threads for the one Prokka run')
    args = parser.parse_args()
    if args.fasta is None:
        args.print_help()
    else:
        run_prokka(args.fasta, args.threads)
        
