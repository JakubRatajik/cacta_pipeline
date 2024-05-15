import random
from typing import List, Tuple, Dict

from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore


with open("../randgen20Mb/20MBgenome.fasta", "r") as handle:
    genome = list(SimpleFastaParser(handle))

with open("../RepBase/repbase_CACTA_viridiplantae.fasta", "r") as handle:
    elements = list(SimpleFastaParser(handle))

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
        tsd = "".join(random.choice(['A', 'C', 'G', 'T']) for _ in range(3))
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


with open("randgen20nonest/RG_withInsertions20nonest.fasta", "w") as handle:
    for (ch_title, ch_seq) in genome:
        handle.write(f">{ch_title}\n{ch_seq}\n")

with open("randgen20nonest/insertedCACTASeqs20nonest.fasta", "w"):
    pass

with open("randgen20nonest/insertedCACTASeqs20nonest.fasta", "a") as handle:
    for (i, (ch_title, ch_seq)) in enumerate(genome):
        handle.writelines(f">{title}, {start}-{end}\n"
                          f"{ch_seq[start:end]}\n"
                          for (title, start, end) in positions[i])

with open("randgen20nonest/insertedCACTASeqs20nonest.gff3", "w") as handle:
    handle.write("##gff-version 3\n")

with open("randgen20nonest/insertedCACTASeqs20nonest.gff3", "a") as handle:
    for (i, (ch_title, ch_seq)) in enumerate(genome):
        for (title, start, end) in positions[i]:
            handle.write(f"{i + 1}\tinsert_cacta.py\tCACTA_TIR_transposon"
                         f"\t{start}\t{end}\t.\t+\t.\tSeqName={title}\n")

print(f"{len(elements)} (x2) CACTA have been inserted into"
      f" genome. {nested_elements_count} of inserted elements are nested.")
