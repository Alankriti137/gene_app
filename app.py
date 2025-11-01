# app.py
# GeneCheck - BRCA1 + 9 other genes (single-file Streamlit app)
# Save as app.py in your gene_app folder. Put images in gene_app/assets/
# Date: Oct 2025 (used as dataset timestamp in the UI/disclaimer)

import streamlit as st
import pandas as pd
import numpy as np
import os
from datetime import datetime
from textwrap import dedent

st.set_page_config(
    page_title="GeneCheck",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto",
)

# ---------- Metadata / Disclaimer ----------
DATA_TIMESTAMP = "October 2025"
DISCLAIMER = (
    "Data compiled from public resources (NCBI, UniProt, RCSB, Human Protein Atlas). "
    "For educational and informational purposes only — not for medical decision-making."
)

st.markdown("<style> .card-title{font-size:34px; font-weight:700;} </style>", unsafe_allow_html=True)

# ---------- Built-in gene database (example curated entries) ----------
# For each gene: name-> dict with function, mutations (list), expression (dict tissue->value),
# protein structure id (RCSB PDB id or filename base), detailed text.
GENES = {
    "BRCA1": {
        "name": "BRCA1",
        "function": "DNA repair and tumor suppression (homologous recombination).",
        "mutations": ["185delAG", "5382insC", "C61G", "R1443X"],
        "expression": {"Breast": 90, "Ovary": 85, "Prostate": 60, "Breast (normal)": 30},
        "protein_image": "BRCA1_structure.png",
        "protein_pdb": "1JNX",
        "summary": dedent("""\
            BRCA1 helps repair double-strand DNA breaks by homologous recombination and
            helps maintain genome stability. Pathogenic variants strongly increase lifetime risk
            of breast and ovarian cancer. Clinical testing and genetic counseling recommended when
            family history is suggestive.
        """),
        "more": "Key reviews: NCBI Gene / OMIM / clinical resources.",
    },
    "TP53": {
        "name": "TP53",
        "function": "Tumor suppressor, guardian of the genome; transcription factor controlling cell cycle and apoptosis.",
        "mutations": ["R175H", "R248Q", "R273H"],
        "expression": {"Colon": 70, "Lung": 80, "Liver": 65},
        "protein_image": "TP53_structure.png",
        "protein_pdb": "1TUP",
        "summary": "TP53 encodes p53, a master regulator of the DNA damage response and apoptosis.",
        "more": "TP53 mutations are common across many cancers and are often somatic.",
    },
    "CFTR": {
        "name": "CFTR",
        "function": "Chloride ion channel involved in fluid transport across epithelia.",
        "mutations": ["ΔF508", "G551D"],
        "expression": {"Lung": 95, "Pancreas": 90},
        "protein_image": "CFTR_structure.png",
        "protein_pdb": "5UAK",
        "summary": "CFTR variants cause cystic fibrosis; ΔF508 is the most frequent pathogenic variant.",
        "more": "Therapies exist targeting specific CFTR mutations (e.g., modulators).",
    },
    "EGFR": {
        "name": "EGFR",
        "function": "Receptor tyrosine kinase that promotes cell growth signaling.",
        "mutations": ["L858R", "T790M"],
        "expression": {"Lung": 85, "Brain": 40, "Skin": 30},
        "protein_image": "EGFR_structure.png",
        "protein_pdb": "2ITY",
        "summary": "EGFR mutations can drive certain lung cancers; targeted inhibitors exist.",
        "more": "EGFR mutation testing is standard for non-small cell lung cancer.",
    },
    "APOE": {
        "name": "APOE",
        "function": "Lipid transport and cholesterol metabolism.",
        "mutations": ["E2", "E3", "E4 variants"],
        "expression": {"Brain": 70, "Liver": 90},
        "protein_image": "APOE_structure.png",
        "protein_pdb": "1GS9",
        "summary": "APOE alleles (especially APOE-e4) affect Alzheimer's disease risk and lipid metabolism.",
        "more": "APOE genotype is a risk factor, not a deterministic predictor.",
    },
    "BRCA2": {
        "name": "BRCA2",
        "function": "DNA repair — homologous recombination (BRCA2 interacts with RAD51).",
        "mutations": ["6174delT", "999del5"],
        "expression": {"Breast": 88, "Ovary": 82, "Prostate": 55},
        "protein_image": "BRCA2_structure.png",
        "protein_pdb": "1MIU",
        "summary": "Like BRCA1, BRCA2 variants increase risk of hereditary breast and ovarian cancer.",
        "more": "",
    },
    "KRAS": {
        "name": "KRAS",
        "function": "Small GTPase involved in signal transduction (MAPK pathway).",
        "mutations": ["G12D", "G13D", "Q61H"],
        "expression": {"Pancreas": 85, "Colon": 80},
        "protein_image": "KRAS_structure.png",
        "protein_pdb": "4LUC",
        "summary": "KRAS activating mutations are common driver mutations in many cancers.",
        "more": "",
    },
    "MYC": {
        "name": "MYC",
        "function": "Transcription factor that drives proliferation and growth.",
        "mutations": ["Amplification", "T58A"],
        "expression": {"Lymph": 75, "Liver": 60},
        "protein_image": "MYC_structure.png",
        "protein_pdb": "6G6K",
        "summary": "MYC is often dysregulated in cancer via amplification or transcriptional upregulation.",
        "more": "",
    },
    "LDLR": {
        "name": "LDLR",
        "function": "Low-density lipoprotein receptor: mediates uptake of cholesterol.",
        "mutations": ["D206E", "C646Y"],
        "expression": {"Liver": 95, "Adrenal": 45},
        "protein_image": "LDLR_structure.png",
        "protein_pdb": "3P5B",
        "summary": "LDLR variants cause familial hypercholesterolemia in many cases.",
        "more": "",
    },
    "MTHFR": {
        "name": "MTHFR",
        "function": "Folate metabolism enzyme (methylene tetrahydrofolate reductase).",
        "mutations": ["C677T", "A1298C"],
        "expression": {"Liver": 70, "Kidney": 50},
        "protein_image": "MTHFR_structure.png",
        "protein_pdb": "6FCX",
        "summary": "MTHFR polymorphisms can affect folate metabolism; clinical significance varies.",
        "more": "",
    },
}

# ---------- Helper functions ----------
ASSETS_DIR = os.path.join(os.getcwd(), "assets")  # expects gene_app/assets/

def asset_exists(filename: str) -> bool:
    if not filename:
        return False
    return os.path.isfile(os.path.join(ASSETS_DIR, filename))

def pdb_link(gene_entry):
    pdb = gene_entry.get("protein_pdb")
    if not pdb:
        return None
    return f"https://www.rcsb.org/structure/{pdb}"

def ncbi_gene_link(gene_symbol):
    return f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_symbol}"

def uniprot_link(gene_symbol):
    return f"https://www.uniprot.org/uniprot/?query={gene_symbol}&sort=score"

def find_gene(symbol: str):
    s = symbol.strip().upper()
    return GENES.get(s)

def make_expression_df(expr_dict):
    # convert tissues -> numeric values (0-100)
    df = pd.DataFrame(list(expr_dict.items()), columns=["Tissue", "Expression"])
    df = df.sort_values("Expression", ascending=False)
    return df

# ---------- Layout ----------
st.title("🧬 GeneCheck")
st.write(f"Explore gene details – function, known mutations, expression levels, and protein structures. (Data snapshot: **{DATA_TIMESTAMP}**)")

st.caption(DISCLAIMER)

# Search area
st.markdown("### 🔎 Search for a gene")
col1, col2 = st.columns([3, 1])
with col1:
    gene_input = st.text_input("Enter gene name or symbol (e.g., BRCA1, TP53, CFTR):", value="BRCA1")
with col2:
    _ = st.button("Search")

if not gene_input:
    st.info("Type a gene symbol above (for example: BRCA1) and press Search.")
    st.stop()

# Fetch gene
entry = find_gene(gene_input)
if not entry:
    # try fuzzy: look for gene symbol substring
    matches = [g for g in GENES.keys() if gene_input.strip().upper() in g]
    if matches:
        entry = GENES[matches[0]]
    else:
        st.error(f"No data found for '{gene_input}'. Try one of: {', '.join(sorted(GENES.keys()))}")
        st.stop()

# ---------- Header card ----------
st.markdown("---")
hcol1, hcol2 = st.columns([2,3])
with hcol1:
    st.image(
        os.path.join(ASSETS_DIR, entry.get("protein_image")) if asset_exists(entry.get("protein_image")) else "https://via.placeholder.com/360x360?text=Protein",
        width=280
    )
with hcol2:
    st.markdown(f"<div class='card-title'>{entry['name']}</div>", unsafe_allow_html=True)
    st.write(f"**Function:** {entry['function']}")
    st.write(f"**Summary:** {entry['summary']}")
    st.write("**Quick links:** ", unsafe_allow_html=True)
    links = []
    links.append(f"[NCBI Gene]({ncbi_gene_link(entry['name'])})")
    links.append(f"[UniProt]({uniprot_link(entry['name'])})")
    if pdb_link(entry):
        links.append(f"[RCSB PDB]({pdb_link(entry)})")
    st.markdown(" • ".join(links))

st.markdown("---")

# ---------- Info columns ----------
c1, c2 = st.columns([2,3])
with c1:
    st.subheader("Gene at a glance")
    st.table(pd.DataFrame({
        "Field": ["Name", "Function", "Protein structure"],
        "Value": [entry["name"], entry["function"], entry.get("protein_pdb", "—")]
    }))

    st.subheader("Known mutations")
    if entry.get("mutations"):
        st.write(", ".join(entry["mutations"]))
    else:
        st.write("No curated mutations in this dataset.")

    st.subheader("Protein structure")
    if asset_exists(entry.get("protein_image")):
        st.image(os.path.join(ASSETS_DIR, entry["protein_image"]), caption=f"{entry['name']} structure (asset)", use_column_width=True)
    elif entry.get("protein_pdb"):
        st.info("3D structure available on RCSB PDB.")
        st.write(f"[Open {entry['protein_pdb']} on RCSB]({pdb_link(entry)})")
    else:
        st.write("No local structure available.")

with c2:
    st.subheader("Expression levels (selected tissues)")
    expr = entry.get("expression", {})
    if expr:
        df_expr = make_expression_df(expr)
        st.bar_chart(data=df_expr.set_index("Tissue"))
        # show table
        st.dataframe(df_expr, width=600, height=160)
    else:
        st.write("No expression data available.")

    st.subheader("Interactive explorer")
    with st.expander("View mutations mapped to example gene regions (illustrative)"):
        # We avoid chromosome-scale plotting; provide an illustrative schematic using text + bar_chart
import altair as alt

muts = entry.get("mutations", [])
if muts:
    df_m = pd.DataFrame({"Mutation": muts, "Count": [1]*len(muts)})
    chart = (
        alt.Chart(df_m)
        .mark_bar()
        .encode(x="Mutation:N", y="Count:Q")
        .properties(width="container")
    )
    st.altair_chart(chart, use_container_width=True)
else:
    st.write("No mutation list to show here.")
            st.write("Note: this is an illustrative marker plot — not genomic coordinates. If you want real chromosome coordinates, I can add those when we fetch genome annotations.")
        else:
            st.write("No mutation list to show here.")

# ---------- Deep-dive / Expanders ----------
st.markdown("---")
st.subheader("More details")
with st.expander("Clinical and research notes"):
    st.write(entry.get("more") or "No extra notes available.")
    st.markdown(
        dedent(f"""
        **Credible resources**  
        - NCBI Gene: {ncbi_gene_link(entry['name'])}  
        - UniProt search: {uniprot_link(entry['name'])}  
        - RCSB (if available): {pdb_link(entry) if pdb_link(entry) else '—'}
        """
        )
    )

with st.expander("Related genes in this dataset"):
    links = []
    for sym in sorted(GENES.keys()):
        links.append(f"- **{sym}**: {GENES[sym]['function']}")
    st.write("\n".join(links))

# ---------- Footer / navigation ----------
st.markdown("---")
cols = st.columns([1,1,6])
with cols[0]:
    st.write(f"Data snapshot: **{DATA_TIMESTAMP}**")
with cols[1]:
    st.write(f"App: GeneCheck")
with cols[2]:
    st.caption(DISCLAIMER)

# ---------- End ----------

