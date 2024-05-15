import argparse


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for tir_information.py script.
    :return: Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="tir_information.py",
        description="Retrieves information about TIR of Subclass I DNA"
                    " transposable elements.")

    parser.add_argument("-i",
                        "--in-file",
                        type=str,
                        help="Input FASTA file with transposable element"
                             "sequences",
                        required=True)

    parser.add_argument("-o",
                        "--out-file",
                        type=str,
                        help="Output FASTA file containing transposable"
                             "element sequences with appended TIR information",
                        required=True)

    parser.add_argument("-t",
                        "--tir-length",
                        type=int,
                        default=28,
                        help="TIR length to be aligned")

    args = parser.parse_args()

    if args.in_file is None or args.out_file is None:
        parser.error("-i and -o parameters must be specified")

    return args
