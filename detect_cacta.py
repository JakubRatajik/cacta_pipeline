import argparse
import tempfile
import time
from typing import List, Tuple, IO
from typing_extensions import TypeAlias

from Bio.Seq import Seq  # type: ignore
from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore

from utils import tir_information
from parsing.detect_cacta import parse_arguments


Candidate: TypeAlias = Tuple[str, str, int, int, int]
Record: TypeAlias = Tuple[str, str]
TirPair: TypeAlias = Tuple[int, int]
PositionHash: TypeAlias = Tuple[int, int]

BASE_TO_NUM = {"A": 0, "C": 1, "G": 2, "T": 3}
candidate_id = 1


def tir_tsd_hash(tir: str, tsd: str, opening_tir: bool) -> int:
    """
    For two potential matching TIRs, instead of comparing whether each base of
    one TIR corresponds to reverse complement of another TIR, number encoding
    TIR and TSD (hash) can be calculated and then for each potential matching
    TIR pair only the corresponding hashes can be compared.

    Nucleotide sequence can be unambiguously hashed in quaternary numeral
    system. Each value at given position represents base at that position.
    A: 0, C: 1, G: 2, T: 3
    For closing pair, reverse complement of TIR and direct TSD is hashed, so
    tir_tsd_hash("CACTA", "ATT", false) == tir_tsd_hash("TAGTG", "ATT", true)

    E.g. opening TIR "TAGCT" and corresponding TSD "TTC" is represented in
    quaternary numeral system as 13331203.

    :param tir: TIR sequence
    :param tsd: TSD sequence
    :param opening_tir: True if opening TIR, False otherwise
    :return: TIR + TSD representation in quaternary numeral system
    """

    hash_code = 0
    exponent = 0

    if not opening_tir:
        tir = Seq(tir).reverse_complement()

    for base in tir:
        hash_code += BASE_TO_NUM[base] * (4 ** exponent)
        exponent += 1

    for base in tsd:
        hash_code += BASE_TO_NUM[base] * (4 ** exponent)
        exponent += 1

    return hash_code


def tir_tsd_condensed(tir_position: int, genome_length: int,
                      opening_tir: bool) -> bool:
    """
    Ensures that TIR and TSD are not too close to start or to end of genome.

    :param tir_position: position of TIR in genome
    :param genome_length: length of genome
    :param opening_tir:
    :return:
    """
    return (opening_tir
            and (tir_position < 3 or genome_length < tir_position + 10)) \
        or (not opening_tir
            and (genome_length < tir_position + 8 or tir_position < 5))


def create_prefix_table(pattern: str) -> List[int]:
    pattern_length = len(pattern)
    prefix_table = [0] * pattern_length
    match_count = 0

    for i in range(1, pattern_length):
        while match_count > 0 and pattern[match_count] != pattern[i]:
            match_count = prefix_table[match_count - 1]

        if pattern[match_count] == pattern[i]:
            match_count += 1
            prefix_table[i] = match_count

    return prefix_table


# TODO: should not include occurrences in the name but rather tir, CACTA..
#  but then it should also not have the pattern param...
#  there is no check if the the fasta file is valid. it requires that
#  sequences contain only A, C, T, G, and N letters
#
def find_all_tirs(tir: str, sequence: str,
                  opening_tir: bool) -> List[PositionHash]:
    """
    Employs Knuth-Morris-Pratt algorithm to find all occurrences of TIR in
    genome.
    The search:
        - expects input sequence to contain exclusively uppercase characters
        - does not support entire IUPAC nucleotide code when looking for TIR
            and TSD, i.e. can only contain: A, C, G, T characters.

    :param tir: TIR to be searched for
    :param sequence: nucleotide sequence to be searched in
    :param opening_tir: True if the TIR is opening (at the 5' end),
        False otherwise
    :return: List containing for each occurrence its position in genome
        along with hash value of TIR remainder with corresponding TSD
    """
    sequence_length = len(sequence)
    tir_length = len(tir)
    occurrences: List[PositionHash] = []

    if tir_length == 0:
        return occurrences

    prefix_table = create_prefix_table(tir)

    j = 0

    for i in range(sequence_length):
        while j > 0 and sequence[i] != tir[j]:
            j = prefix_table[j - 1]

        if sequence[i] != tir[j]:
            continue

        j += 1

        if j != tir_length:
            continue

        j = prefix_table[j - 1]
        tir_position = i - tir_length + 1

        if tir_tsd_condensed(tir_position, sequence_length, opening_tir):
            continue

        if opening_tir:
            tir_remainder = sequence[tir_position + 5: tir_position + 10]
            tsd = sequence[tir_position - 3: tir_position]
        else:
            tir_remainder = sequence[tir_position - 5: tir_position]
            tsd = sequence[tir_position + 5: tir_position + 8]

        # Do TIR remainder or TSD contain restricted chars?
        if set(tir_remainder + tsd) - {"A", "C", "G", "T"} != set():
            continue

        hash_code = tir_tsd_hash(tir_remainder, tsd, opening_tir)

        occurrences.append((tir_position, hash_code))

    return occurrences


# For opening_tirs[i] at position X in genome, there are N closing
# TIRs at lower position than X. For opening_tirs[i+1] we can certainly
# skip at least first N closing TIRs.
def filter_matching_tirs(opening_tirs: List[PositionHash],
                         closing_tirs: List[PositionHash],
                         min_length: int, max_length: int) -> List[TirPair]:
    """
    Finds all matching TIR pairs (opening with closing) based on the TIR and
    TSD sequences, as well as their relative position in genome.
    :param opening_tirs: list of position opening TIRs (position and hash)
    :param closing_tirs: list of position closing TIRs (position and hash)
    :param min_length: minimum TE length
    :param max_length: maximum TE length
    :return: list of closing and opening TIR pair positions
    """
    smallest_relevant_closing_tir = 0
    matching_tirs = []
    closing_tirs_len = len(closing_tirs)

    for opening_tir in opening_tirs:
        opening_position, opening_hash = opening_tir[0], opening_tir[1]
        for i in range(smallest_relevant_closing_tir, closing_tirs_len):
            closing_position, closing_hash = closing_tirs[i][0], \
                closing_tirs[i][1]

            if opening_position + max_length < closing_position:
                break

            if closing_position + 4 < opening_position + min_length:
                smallest_relevant_closing_tir += 1
                continue

            if opening_hash == closing_hash:
                matching_tirs.append((opening_position, closing_position))

    return matching_tirs


def retrieve_candidates(record: Record, matching_tirs: List[TirPair],
                        tir_info: bool, elements: List[Candidate],
                        seq_id: int) -> None:
    global candidate_id

    record_title = record[0]
    record_seq = record[1]

    for (opening, closing) in matching_tirs:
        title = f"{record_title}_CACTA{candidate_id}"
        sequence = record_seq[opening:closing + 5]

        if tir_info:
            title, sequence = tir_information.extract_tir_info(
                (title, sequence), tir_information.default_aligner())

        elements.append((title, sequence, opening, closing + 5, seq_id))
        candidate_id += 1


def detect_cacta(record: Record, args: argparse.Namespace,
                 elements: List[Candidate], seq_id: int) -> None:
    print(f"Getting opening CACTA TIRs.")
    cacta_tirs = find_all_tirs("CACTA", record[1], True)

    print(f"Getting closing TAGTG TIRs.")
    tagtg_tirs = find_all_tirs("TAGTG", record[1], False)

    matching_tirs = filter_matching_tirs(cacta_tirs, tagtg_tirs,
                                         args.min_len, args.max_len)

    retrieve_candidates(record, matching_tirs, args.tir_info, elements, seq_id)


def detect_cactg(record: Record, args: argparse.Namespace,
                 elements: List[Candidate], seq_id: int) -> None:
    print(f"Getting opening CACTG TIRs.")
    cactg_tirs = find_all_tirs("CACTG", record[1], True)

    print(f"Getting closing CAGTG TIRs.")
    cagtg_tirs = find_all_tirs("CAGTG", record[1], False)

    matching_tirs = filter_matching_tirs(cactg_tirs, cagtg_tirs,
                                         args.min_len, args.max_len)

    retrieve_candidates(record, matching_tirs, args.tir_info, elements, seq_id)


def detect_all_candidates(temp_handle: IO[str],
                          args: argparse.Namespace) -> List[Candidate]:
    """
    Detects CACTA TE candidates (both CACTA and CACTG) in input DNA sequence.
    :param temp_handle: temporary file handle
    :param args: parsed command line arguments
    :return: list of detected CACTA TE candidates
    """
    candidates: List[Candidate] = []

    for seq_id, record in enumerate(SimpleFastaParser(temp_handle), 1):
        print(f"Processing '{record[0]}' sequence.")

        count_before = len(candidates)

        detect_cacta(record, args, candidates, seq_id)
        cacta_count = len(candidates) - count_before

        detect_cactg(record, args, candidates, seq_id)
        cactg_count = len(candidates) - count_before - cacta_count

        print(f"Found {cacta_count} CACTA-TAGTG and "
              f"{cactg_count} CACTG-CAGTG matching TIRs in {record[0]}.\n")

    print(f"\nOverall, {candidate_id - 1} CACTA candidates were detected.\n")

    return candidates


def file_to_uppercase_temp(temp_handle: IO[str], in_file_name: str) -> None:
    """
    Converts input file to uppercase and writes it to temporary file.
    :param temp_handle: temporary file handle
    :param in_file_name: name of input file
    """
    with open(in_file_name, "r") as in_file_handle:
        for record in SimpleFastaParser(in_file_handle):
            temp_handle.write(f">{record[0]}\n{record[1].upper()}\n")
    temp_handle.seek(0)


def export_fasta(candidates: List[Candidate], out_file: str) -> None:
    print("\nExporting transposon sequences in FASTA format.")

    with open(out_file, "w") as handle:
        for candidate in candidates:
            handle.write(f">{candidate[0]}\n{candidate[1]}\n")

    print(f"Transposon sequences are stored in '{out_file}'")


def export_gff3(candidates: List[Candidate], out_file: str) -> None:
    print("\nExporting transposon annotation in GFF3 format.")

    with open(out_file, "w") as handle:
        handle.write("##gff-version 3\n")

        for title, seq, start, end, seq_id in candidates:
            handle.write(f"{seq_id}\t"
                         f"detect_cacta.py\t"
                         f"CACTA_TIR_transposon\t"
                         f"{start}\t"
                         f"{end}\t"
                         f".\t"
                         f"+\t"
                         f".\t"
                         f"SeqName={title};\n")

    print(f"Transposon sequences are stored in '{out_file}'")


def main() -> None:
    start_time = time.time()
    args = parse_arguments()

    with tempfile.TemporaryFile(mode="r+") as temp_handle:
        file_to_uppercase_temp(temp_handle, args.in_file)
        candidates = detect_all_candidates(temp_handle, args)

    if args.fasta_out:
        export_fasta(candidates, args.fasta_out)
    if args.gff3:
        export_gff3(candidates, args.gff3)

    elapsed_time = time.gmtime(time.time() - start_time)
    print(f"\nDetect CACTA has just finished. Elapsed time: "
          f"{time.strftime('%Hh:%Mm:%Ss', elapsed_time)}.\n")


if __name__ == '__main__':
    main()
