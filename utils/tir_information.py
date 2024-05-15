from typing import Tuple

from Bio import Align  # type: ignore
from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore

from parsing.tir_information import parse_arguments


def extract_tir_info(record: Tuple[str, str], aligner: Align.PairwiseAligner,
                     length: int = 28) -> Tuple[str, str]:
    """
    Retrieve TIR information for trasnposon sequence.
    :param record: tuple with transposon sequence title and sequence
    :param aligner: pairwise aligner object
    :param length: TIR length to be aligned
    """
    title = record[0]
    sequence = record[1]

    best_alignment = max(aligner.align(str.upper(sequence[0:length]),
                                       str.upper(sequence[-length:]),
                                       strand="-"))

    gaps, identities, mismatches = best_alignment.counts()

    new_title = (f"{title}_{gaps + identities + mismatches}"
                 f"bpTIR(m={mismatches}, g={gaps})")

    return new_title, sequence


def default_aligner() -> Align.PairwiseAligner:
    """
    Default aligner with match=2, mismatch=-3, open=-5, extend=-2
    """
    return init_aligner(2, -3, -5, -2)


def init_aligner(match_score: int, mismatch_score: int, open_gap_score: int,
                 extend_gap_score: int) -> Align.PairwiseAligner:
    """
    Set alignment scores for the aligner object
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    return aligner


def main() -> None:
    args = parse_arguments()
    in_file = args.in_file
    out_file = args.out_file
    length = args.tir_length

    with open(in_file, "r") as handle:
        records = list(SimpleFastaParser(handle))

    aligner = default_aligner()

    with_tir_info = map(lambda rec: extract_tir_info(rec, aligner, length),
                        records)

    with open(out_file, "w") as handle:
        for record in with_tir_info:
            handle.write(f">{record[0]}\n{record[1]}\n")


if __name__ == '__main__':
    main()
