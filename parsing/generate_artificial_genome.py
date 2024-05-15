import argparse

from parsing import parsing_utils

"""
This module contains functions for parsing command line
arguments for the generate_artificial_genome.py script.
"""


def validate_genome_size(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 0, 10000000000)


def validate_chromosome_number(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 0, 300)


def validate_chunk_size(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 20, 100000000)


def validate_gc_content(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 0, 100)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for generate_artificial_genome.py script.
    :return: namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="generate_artificial_genome.py",
        description="Tool for generating pseudo-random genomic sequence.")

    parser.add_argument("-o",
                        "--out-file",
                        type=str,
                        help="File to write FASTA formatted sequence into",
                        required=True)

    parser.add_argument("-s",
                        "--size",
                        type=validate_genome_size,
                        help="Size of the genomic sequence to generate (bp)",
                        required=True)

    parser.add_argument("-n",
                        "--num-of-chromosomes",
                        type=validate_chromosome_number,
                        help="Number of genome chromosomes.",
                        required=True)

    parser.add_argument("-gc",
                        "--gc-content",
                        type=validate_gc_content,
                        help="GC content (%)",
                        required=True)

    parser.add_argument("-ch",
                        "--chunk-size",
                        type=validate_chunk_size,
                        help="Iteratively write genomic sequence of chunk size"
                             " to file not overfilling memory",
                        required=True)

    return parser.parse_args()
