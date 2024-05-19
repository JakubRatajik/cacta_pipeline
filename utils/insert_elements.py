from parsing.insert_elements import parse_arguments
import random
from typing import List, Tuple, Dict
from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore


def read_genome(file_path: str) -> List[Tuple[str, str]]:
    """
    Load genome from file.
    :param file_path: path to file with genome
    :return: list of chromosome records
    """
    with open(file_path, "r") as handle:
        return list(SimpleFastaParser(handle))


def read_elements(file_path: str) -> List[Tuple[str, str]]:
    """
    Load transposable elements from file.
    :param file_path: path to file with transposable elements
    :return: list of transposable element records
    """
    with open(file_path, "r") as handle:
        return list(SimpleFastaParser(handle))


def write_genome(file_path: str, genome: List[Tuple[str, str]]) -> None:
    """
    Write genome with inserted elements to file.
    :param file_path: path to output file
    :param genome: list of chromosome records
    """
    with open(file_path, "w") as handle:
        for (ch_title, ch_seq) in genome:
            handle.write(f">{ch_title}\n{ch_seq}\n")


def write_inserted_sequences(file_path: str, genome: List[Tuple[str, str]],
                             positions: Dict[
                                 int, List[Tuple[str, int, int]]]) -> None:
    """
    Write inserted elements to file.
    :param file_path: file to write sequences to
    :param genome: list of chromosome records
    :param positions: positions of inserted elements
    """
    with open(file_path, "w"):
        pass

    with open(file_path, "a") as handle:
        for (i, (ch_title, ch_seq)) in enumerate(genome):
            handle.writelines(f">{title}, {start}-{end}\n"
                              f"{ch_seq[start:end]}\n"
                              for (title, start, end) in positions[i])


def write_gff3(file_path: str, genome: List[Tuple[str, str]],
               positions: Dict[int, List[Tuple[str, int, int]]]) -> None:
    """
    Write annotation of inserted elements.
    :param file_path: path to output file
    :param genome: list of chromosome records
    :param positions: positions of inserted elements
    """
    with open(file_path, "w") as handle:
        handle.write("##gff-version 3\n")

    with open(file_path, "a") as handle:
        for (i, _) in enumerate(genome):
            for (title, start, end) in positions[i]:
                handle.write(f"{i + 1}\tinsert_cacta.py\tCACTA_TIR_transposon"
                             f"\t{start}\t{end}\t.\t+\t.\tSeqName={title}\n")


def insert_elements(genome: List[Tuple[str, str]],
                    elements: List[Tuple[str, str]], output_dir: str) -> None:
    """
    Insert transposable elements into genome.
    :param genome: list of chromosome records
    :param elements: transposable elements to be inserted
    :param output_dir: output directory
    :return:
    """
    nested_elements_count = 0
    positions: Dict[int, List[Tuple[str, int, int]]] = \
        {i: [] for i in range(len(genome))}

    for element in elements:
        element_name, element_seq = element
        element_length = len(element_seq)

        for j in range(1, 3):
            chromosome = random.randrange(0, len(genome))
            ch_title, ch_seq = genome[chromosome]
            insert_position = random.randrange(0, len(ch_seq))
            tsd = "".join(
                random.choice(['A', 'C', 'G', 'T']) for _ in range(3))
            element_start = insert_position + len(tsd)
            element_end = element_start + element_length
            insertion_length = element_length + len(tsd) * 2

            genome[chromosome] = (
                ch_title, (ch_seq[:insert_position] + tsd + element_seq
                           + tsd + ch_seq[insert_position:]))

            ch_positions = positions[chromosome]
            for i in range(len(ch_positions)):
                (title, start, end) = ch_positions[i]

                if end <= insert_position:
                    ch_positions[i] = (title, start, end)
                elif insert_position <= start:
                    ch_positions[i] = (title, start + insertion_length,
                                       end + insertion_length)
                else:
                    ch_positions[i] = (title, start, end + insertion_length)
                    nested_elements_count += 1

            ch_positions.append((f"{element[0]}_{j}_ch{chromosome + 1}",
                                 element_start, element_end))

    write_genome(f"{output_dir}/genome_with_insertions.fasta", genome)
    write_inserted_sequences(f"{output_dir}/inserted_elements.fasta",
                             genome, positions)
    write_gff3(f"{output_dir}/inserted_elements.gff3", genome,
               positions)

    print(f"{len(elements)} (x2) CACTA have been inserted into"
          f" genome. {nested_elements_count} of inserted elements are nested.")


def main() -> None:
    args = parse_arguments()

    genome = read_genome(args.genome)
    elements = read_elements(args.elements)
    insert_elements(genome, elements, args.output_dir)


if __name__ == '__main__':
    main()
