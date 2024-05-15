import numpy as np

from parsing.generate_artificial_genome import parse_arguments


def generate_artificial_sequence(size: int, chromosomes: int, gc_content: int,
                                 chunk_size: int, out_file: str) -> None:
    """
    Generate artificial genomic sequence with specified size, number of
    chromosomes and GC content.
    :param size: overall size of the genome
    :param chromosomes: number of chromosomes
    :param gc_content: GC content of the genome
    :param chunk_size: chunk size to write to file iteratively
    :param out_file: file to write the sequence into
    """
    nts = ["A", "C", "G", "T"]
    at_content = 100 - gc_content
    a_probability = t_probability = at_content / 200
    c_probability = g_probability = gc_content / 200
    probabilities = [a_probability, c_probability,
                     g_probability, t_probability]

    size = (size // chromosomes) * chromosomes

    print(f"Generating {size}bp genomic sequence within {chromosomes}"
          f" chromosomes and {gc_content}% GC content")

    if size < 1000000:
        str_size = f"{size}bp"
    elif size < 1000000000:
        str_size = f"{round(size / 1000000, 2)}Mb"
    else:
        str_size = f"{round(size / 1000000000, 2)}Gb"

    with open(out_file, "w"):
        pass

    for i in range(1, chromosomes + 1):
        chromosome_size = size // chromosomes
        print(f"Generating chromosome {i} of {chromosomes}"
              f" ({chromosome_size}bp)")

        with open(out_file, "a") as handle:
            handle.write(f">randomGenome_{str_size}_{gc_content}GC_ch{i}\n")

        while chromosome_size > 0:
            ch = min(chromosome_size, chunk_size)
            chromosome_size -= ch

            sequence = "".join(np.random.choice(nts, size=ch, p=probabilities))

            with open(out_file, "a") as handle:
                handle.write(f"{sequence}")

        with open(out_file, "a") as handle:
            handle.write("\n")

    print(f"Genomic sequence has been written to '{out_file}'")


def main() -> None:
    args = parse_arguments()

    generate_artificial_sequence(args.size, args.num_of_chromosomes,
                                 args.gc_content, args.chunk_size,
                                 args.out_file)


if __name__ == '__main__':
    main()
