import argparse


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for the insert_elements.py script.
    :return: Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="insert_elements.py",
        description="Insert transposable elements into a genome.")

    parser.add_argument("-g",
                        "--genome",
                        type=str,
                        help="Path to the FASTA file with genome sequence"
                             "to insert transposons into.",
                        required=True)

    parser.add_argument("-e",
                        "--elements",
                        type=str,
                        help="Path to the FASTA file containing the"
                             "transposon sequences.",
                        required=True)

    parser.add_argument("-o",
                        "--output-dir",
                        type=str,
                        help="The output directory path.",
                        required=True)

    return parser.parse_args()
