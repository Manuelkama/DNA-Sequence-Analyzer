from Bio import SeqIO
from Bio.Data import CodonTable


def analyze_sequence(sequence):
    sequence = "".join(str(sequence).split()).upper()

    bases = ["A", "T", "C", "G"]
    counts = {base: sequence.count(base) for base in bases}

    total_length = len(sequence)
    canonical_total = sum(counts.values())

    gc_content = (
        (counts["G"] + counts["C"]) / canonical_total * 100
        if canonical_total else 0
    )

    at_content = (
        (counts["A"] + counts["T"]) / canonical_total * 100
        if canonical_total else 0
    )

    valid_iupac = set("ATCGRYSWKMBDHVN")

    ambiguous = sum(
        1 for base in sequence
        if base in valid_iupac and base not in bases
    )

    invalid = sum(
        1 for base in sequence
        if base not in valid_iupac
    )

    return {
        "Sequence length": total_length,
        "Base counts": counts,
        "GC content (%)": round(gc_content, 2),
        "AT content (%)": round(at_content, 2),
        "Ambiguous bases": ambiguous,
        "Invalid characters": invalid
    }


def extract_genes(record):
    genes = {}

    for feature in record.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get(
                "gene", ["Unknown"]
            )[0]

            sequence = feature.extract(record.seq)

            genes[gene] = {
                "product": feature.qualifiers.get(
                    "product", ["Unknown"]
                )[0],
                "sequence": sequence,
                "length": len(sequence),
                "strand": feature.location.strand,
                "transl_table": int(
                    feature.qualifiers.get(
                        "transl_table", ["1"]
                    )[0]
                )
            }

    return genes


def translate_gene(gene_data):
    sequence = gene_data["sequence"]
    table_id = gene_data["transl_table"]

    complete_length = (len(sequence) // 3) * 3
    coding_sequence = sequence[:complete_length]

    protein = str(
        coding_sequence.translate(table=table_id)
    )

    genetic_code = CodonTable.unambiguous_dna_by_id[table_id]

    if str(coding_sequence[:3]) in genetic_code.start_codons:
        protein = "M" + protein[1:]

    if protein.endswith("*"):
        protein = protein[:-1]

    return protein


def analyze_genome(fasta_path, genbank_path=None):
    record = SeqIO.read(str(fasta_path), "fasta")

    results = {
        "sequence_id": record.id,
        "statistics": analyze_sequence(record.seq),
        "genes": {}
    }

    if genbank_path:
        annotated_record = SeqIO.read(
            str(genbank_path), "genbank"
        )

        genes = extract_genes(annotated_record)

        for gene, data in genes.items():
            protein = translate_gene(data)

            results["genes"][gene] = {
                "product": data["product"],
                "length": data["length"],
                "strand": data["strand"],
                "protein_length": len(protein),
                "protein_sequence": protein
            }

    return results
