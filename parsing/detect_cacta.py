import argparse

from parsing import parsing_utils


"""
This module contains functions for parsing command line
arguments for the detect_cacta.py script.
"""


def validate_min_length(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 50, 23018)


def validate_max_length(arg: str) -> int:
    return parsing_utils.validate_arg_bounds(arg, 50, 30000)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for detect_cacta.py script.
    :return: Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="detect_cacta.py",
        description="A structure-based tool for detecting CACTA transposable"
                    "elements in genomic sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",
                        "--in-file",
                        type=str,
                        help="Input DNA sequence file",
                        required=True)

    parser.add_argument("-fa",
                        "--fasta-out",
                        type=str,
                        help="Output FASTA file with candidate sequences")

    parser.add_argument("-g",
                        "--gff3",
                        type=str,
                        help="Output GFF3 file with candidate annotation")

    parser.add_argument("--min-len",
                        type=validate_min_length,
                        default=50,
                        help="Miminum transposon length")

    parser.add_argument("--max-len",
                        type=validate_max_length,
                        default=23018,
                        help="Maximum transposon length")

    parser.add_argument('--tir-info',
                        action='store_true',
                        help="Append TIR info (TIR length, mismatch count, "
                             "gap count) to candidate name")
    parser.set_defaults(tir_info=False)

    args = parser.parse_args()

    if args.fasta_out is None and args.gff3 is None:
        parser.error("-fa or -gff3 must be specified")

    return args
