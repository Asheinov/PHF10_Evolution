def read_fasta(path):

    sequence = ""

    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()

    return sequence