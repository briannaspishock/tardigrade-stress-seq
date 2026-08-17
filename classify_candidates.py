candidates = {
    "GBZR01011776.1 (Hduj_v1_tr_11787)": {
        "Signal_Peptide": False,
        "Fold_Change": "+18.49",
        "pI": 8.95,
        "Family": "CAHS family (Cytosolic Vitrification Protein)"
    },
    "GBZR01007799.1 (Hduj_v1_tr_07806)": {
        "Signal_Peptide": True,
        "Fold_Change": "+16.70",
        "pI": 10.31,
        "Family": "SAHS candidate (Secretory Stress Protein)"
    },
    "GBZR01000616.1 (Hduj_v1_tr_00616)": {
        "Signal_Peptide": False,
        "Fold_Change": "+13.52",
        "pI": 9.58,
        "Family": "Tardigrade Disordered Protein (TDP / Chaperone)"
    },
    "GBZR01011008.1 (Hduj_v1_tr_11015)": {
        "Signal_Peptide": False,
        "Fold_Change": "+4.43",
        "pI": 5.82,
        "Family": "Repetitive Disordered Linker Protein"
    }
}

print("=" * 74)
print(">> FUNCTIONAL FAMILY CLASSIFICATION (CAHS vs. SAHS)")
print("=" * 74)

for gene, meta in candidates.items():
    sig_status = "Present (Secreted)" if meta.get("Signal_Peptide") else "Absent (Cytosolic)"
    fc = meta.get("Fold_Change")
    pi_val = meta.get("pI")
    fam = meta.get("Family")
    
    print(f"\n[*] {gene}")
    print(f"    - Expression Induction : {fc} Log2FC")
    print(f"    - Predicted Isoelectric: {pi_val}")
    print(f"    - Secretory Signal     : {sig_status}")
    print(f"    - Inferred Family      : {fam}")

print("=" * 74)
