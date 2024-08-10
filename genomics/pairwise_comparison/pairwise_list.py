#!/usr/bin/env python

import sys
import argparse

def usage():
    """Print usage instructions for the script."""
    print("Usage: python pairwise_list.py --entries <entry1> <entry2> ... <entryN>")
    print("Generate all possible pairs from the list of entries.")

def main(argv):
    """Main function to parse arguments and generate pairs."""
    # Define command line options
    CLI = argparse.ArgumentParser(description="Generate all possible pairs from a list of entries.")

    CLI.add_argument(
        "--entries",  # Command line argument name
        nargs="*",  # Accepts zero or more values
        type=str,
        default=[1],  # Default value if no entries are provided
        help="List of entries to pair."
    )

    # Parse the command line arguments
    args = CLI.parse_args()

    # Handle the case when no entries are provided
    if not args.entries:
        usage()
        return

    # Print the list of entries
    print()
    print("Entries: %r" % args.entries)

    # Generate all possible pairs from the list of entries
    res = [(a, b) for idx, a in enumerate(args.entries) for b in args.entries[idx + 1:]]

    # Print the result
    print()
    print("All possible pairs: " + str(res))
    print()

if __name__ == "__main__":
    main(sys.argv[1:])

