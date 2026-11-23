from analysis import (
    analyze_sequence,
    extract_genes,
    translate_gene
)
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import pandas as pd
from io import StringIO

# Page configuration
st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

# Application header
st.title("🧬 DNA Sequence Analyzer")

st.write(
    "Analyze DNA sequences, explore annotated genes, "
    "and translate coding sequences into proteins."
)

st.divider()

# Sequence upload
st.subheader("Upload your sequences")

fasta_file = st.file_uploader(
    "Upload DNA sequence (FASTA)",
    type=["fasta", "fa", "fna"]
)

genbank_file = st.file_uploader(
    "Upload annotated sequence (GenBank)",
    type=["gb", "gbk", "genbank"]
)

if fasta_file:
    st.success(f"FASTA file uploaded: {fasta_file.name}")

if genbank_file:
    st.success(f"GenBank file uploaded: {genbank_file.name}")

if fasta_file:
    fasta_record = SeqIO.read(fasta_file, "fasta")

    results = analyze_sequence(fasta_record.seq)

    st.divider()
    st.header("Sequence Analysis")

    col1, col2, col3 = st.columns(3)

    col1.metric(
        "Sequence Length",
        f"{results['Sequence length']:,} bp"
    )

    col2.metric(
        "GC Content",
        f"{results['GC content (%)']}%"
    )

    col3.metric(
        "AT Content",
        f"{results['AT content (%)']}%"
    )

    st.subheader("Nucleotide Composition")

    st.bar_chart(
        pd.Series(results["Base counts"], name="Count")
    )

    st.subheader("Sequence Quality")

    st.write("Ambiguous bases:", results["Ambiguous bases"])
    st.write("Invalid characters:", results["Invalid characters"])

    if genbank_file:
        genbank_file.seek(0)
        genbank_record = SeqIO.read(genbank_file, "genbank")

        genes = extract_genes(genbank_record)

        st.divider()
        st.header("Gene Analysis")

        selected_gene = st.selectbox(
            "Select a gene",
            list(genes.keys())
        )

        gene_data = genes[selected_gene]
        protein = translate_gene(gene_data)

        st.write("**Gene:**", selected_gene)
        st.write("**Product:**", gene_data["product"])
        st.write("**DNA length:**", gene_data["length"], "bp")
        st.write("**Protein length:**", len(protein), "amino acids")

        with st.expander("View protein sequence"):
            st.code(protein)
