import streamlit as st
import pandas as pd
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from bs4 import BeautifulSoup

# -----------------------------
# PAGE CONFIG
# -----------------------------
st.set_page_config(page_title="Protein Analysis Suite", layout="wide")

# -----------------------------
# DARK UI
# -----------------------------
st.markdown("""
<style>
.stApp {
    background-color: #0e1117;
    color: white;
}

.card {
    background-color: #161b22;
    padding: 20px;
    border-radius: 12px;
    margin-bottom: 15px;
    border: 1px solid #30363d;
}

h1, h2, h3 {
    color: #58a6ff;
}

.team-container {
    display: flex;
    justify-content: center;
    align-items: stretch;
    gap: 40px;
    margin-top: 30px;
    flex-wrap: wrap;
}

.team-card {
    background-color: #161b22;
    padding: 20px 30px;
    border-radius: 12px;
    text-align: center;
    box-shadow: 0 4px 15px rgba(0,0,0,0.25);
    width: 260px;
    border: 1px solid #30363d;
}

.team-card h3 {
    margin-bottom: 10px;
    color: #58a6ff !important;
    font-size: 20px;
}

.team-card p, .team-card a {
    margin: 0;
    font-size: 14px;
    color: white !important;
    text-decoration: none;
}

.team-card a:hover {
    color: #58a6ff !important;
    text-decoration: underline;
}
section[data-testid="stSidebar"] .stRadio label {
    font-size: 22px !important;
}
</style>
""", unsafe_allow_html=True)

# -----------------------------
# CLEAN SEQUENCE
# -----------------------------
def clean_sequence(seq):
    seq = seq.replace(" ", "").replace("\n", "").upper()
    seq = seq.replace("U", "C")   # Selenocysteine -> C
    seq = seq.replace("O", "K")   # Pyrrolysine -> K
    return seq


# -----------------------------
# FETCH SEQUENCE
# -----------------------------
def fetch_sequence(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    r = requests.get(url, timeout=15)

    if r.status_code == 200:
        return "".join(r.text.split("\n")[1:])

    return None


# -----------------------------
# FETCH ISOFORMS
# -----------------------------
def fetch_isoforms(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    r = requests.get(url, timeout=15)

    isoforms = []

    if r.status_code == 200:
        data = r.json()

        for c in data.get("comments", []):
            if c.get("commentType") == "ALTERNATIVE PRODUCTS":

                for iso in c.get("isoforms", []):
                    iso_id = iso.get("isoformIds", [None])[0]

                    if iso_id:
                        fasta_url = f"https://rest.uniprot.org/uniprotkb/{iso_id}.fasta"
                        res = requests.get(fasta_url, timeout=15)

                        if res.status_code == 200:
                            seq = "".join(res.text.split("\n")[1:])
                            isoforms.append((iso_id, seq))

    return isoforms


# -----------------------------
# ANALYSIS
# -----------------------------
def analyze(seq):
    seq = clean_sequence(seq)
    p = ProteinAnalysis(seq)

    return {
        "Molecular Weight": round(p.molecular_weight(), 2),
        "Isoelectric Point": round(p.isoelectric_point(), 2),
        "Instability Index": round(p.instability_index(), 2),
        "GRAVY": round(p.gravy(), 2)
    }


# -----------------------------
# INTERPRETATION
# -----------------------------
def functional_interpretation(res):
    insights = []

    if res["Instability Index"] < 40:
        insights.append("Stable protein")
    else:
        insights.append("Unstable protein")

    if res["GRAVY"] > 0:
        insights.append("Hydrophobic")
    else:
        insights.append("Hydrophilic")

    if res["Isoelectric Point"] > 7:
        insights.append("Basic")
    else:
        insights.append("Acidic")

    return insights


# -----------------------------
# COMPARE ISOFORMS
# -----------------------------
def compare_isoforms(isoforms):
    data = []

    for name, seq in isoforms:
        res = analyze(seq)
        res["Isoform"] = name
        data.append(res)

    return pd.DataFrame(data).set_index("Isoform")


# -----------------------------
# REPORT
# -----------------------------
def validate_sequence(seq):
    lines = seq.strip().splitlines()

    has_header = len(lines) > 0 and lines[0].startswith(">")

    sequence = "".join(
        line.strip()
        for line in lines
        if not line.startswith(">")
    ).upper()

    valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = sorted(set(ch for ch in sequence if ch not in valid_chars))

    return invalid, has_header


def submit_to_server(url, payload):
    r = requests.post(url, data=payload, timeout=60)
    r.raise_for_status()
    return r.text
def generate_report(res, insights):
    text = "Protein Characterization Report\n\n"

    for k, v in res.items():
        text += f"{k}: {v}\n"

    text += "\nInsights:\n"

    for i in insights:
        text += f"- {i}\n"

    return text
# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.title("Navigation")
section = st.sidebar.radio("Go to", ["Tool", "About", "Team"])


# ==================================================
# TOOL PAGE
# ==================================================
if section == "Tool":

    tool = st.sidebar.selectbox(
        "Select Tool",
        ["Protein Analysis", "Isoform Analysis", "SubCellular Location"]
    )

    st.title("Protein Characterization Platform")

    # -------------------------
    # PROTEIN ANALYSIS
    # -------------------------
    if tool == "Protein Analysis":

        mode = st.radio("Input Type", ["Sequence", "UniProt ID"])

        if mode == "Sequence":
            seq = st.text_area("Enter Protein Sequence")
        else:
            uid = st.text_input("Enter UniProt ID")
            seq = fetch_sequence(uid) if uid else None

        if st.button("Analyze Protein"):

            if not seq:
                st.error("Invalid input.")
            else:
                try:
                    res = analyze(seq)
                    insights = functional_interpretation(res)

                    st.markdown('<div class="card">', unsafe_allow_html=True)
                    st.subheader("ProtParam Results")
                    st.dataframe(
                        pd.DataFrame(
                            res.items(),
                            columns=["Property", "Value"]
                        ),
                        use_container_width=True
                    )
                    st.markdown('</div>', unsafe_allow_html=True)

                    st.markdown('<div class="card">', unsafe_allow_html=True)
                    st.subheader("Functional Insights")

                    for i in insights:
                        st.write(f"• {i}")

                    st.markdown('</div>', unsafe_allow_html=True)

                    st.subheader("🔗 External Validation (ProtParam)")
                    st.markdown(
                        "[Open ExPASy ProtParam]"
                        "(https://web.expasy.org/protparam/)"
                    )
                    st.code(clean_sequence(seq), language="text")
                    st.caption(
                        "Copy the sequence above and paste it into "
                        "ExPASy ProtParam for validation."
                    )

                    report = generate_report(res, insights)

                    st.download_button(
                        "📄 Download Report",
                        report,
                        file_name="protein_report.txt"
                    )

                except Exception as e:
                    st.error(f"Analysis failed: {e}")
    # -------------------------
    # ISOFORM ANALYSIS
    # -------------------------
    elif tool == "Isoform Analysis":

        uid = st.text_input("Enter UniProt ID")

        if st.button("Compare Isoforms"):

            if not uid:
                st.warning("Please enter a UniProt ID.")
            else:
                isoforms = fetch_isoforms(uid)

                if len(isoforms) < 2:
                    st.warning("Not enough isoforms found.")
                else:
                    df = compare_isoforms(isoforms)

                    st.subheader("Comparison Table")
                    st.dataframe(df, use_container_width=True)

                    st.subheader("📊 Property-wise Comparison")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.write("Molecular Weight")
                        st.bar_chart(df["Molecular Weight"])

                        st.write("Isoelectric Point (pI)")
                        st.bar_chart(df["Isoelectric Point"])

                    with col2:
                        st.write("Instability Index")
                        st.bar_chart(df["Instability Index"])

                        st.write("GRAVY")
                        st.bar_chart(df["GRAVY"])

                    st.subheader("Key Insights")
                    st.write(f"• Most stable: {df['Instability Index'].idxmin()}")
                    st.write(f"• Most hydrophobic: {df['GRAVY'].idxmax()}")
                    st.write(f"• Largest protein: {df['Molecular Weight'].idxmax()}")
                    st.write(f"• Highest pI: {df['Isoelectric Point'].idxmax()}")

    # -------------------------
    # SUBCELLULAR LOCATION
    # -------------------------
    elif tool == "SubCellular Location":

        st.subheader("Phobius: Topology and Signal Prediction")

        st.radio(
            "Format Option",
            ["Long with Graphics"],
            index=0,
            disabled=True
        )

        seq_phob = st.text_area(
            "Enter Protein Sequence (FASTA with header):",
            key="phob_in"
        )

        if st.button("SUBMIT", key="btn_phob"):

            invalid, has_header = validate_sequence(seq_phob)

            if not has_header:
                st.error("Header required. Example: >protein1")

            elif invalid:
                st.warning(
                    f"Non-amino acid characters detected: {', '.join(invalid)}"
                )

            else:
                try:
                    res = submit_to_server(
                        "https://phobius.sbc.su.se/cgi-bin/predict.pl",
                        {
                            "protseq": seq_phob,
                            "format": "long"
                        }
                    )

                    soup = BeautifulSoup(res, "html.parser")

                    st.write("Prediction Results")

                    pre = soup.find("pre")

                    if pre:
                        text = pre.get_text()
                    else:
                        text = soup.get_text()

                    lines = text.splitlines()
                    rows = []

                    for line in lines:

                        line = line.strip()

                        if line.startswith("FT"):
                            parts = line.split()

                            if len(parts) >= 5:
                                rows.append({
                                    "Feature Type": parts[1],
                                    "Start Position": parts[2],
                                    "End Position": parts[3],
                                    "Subcellular Location": " ".join(parts[4:])
                                })

                    if rows:
                        df = pd.DataFrame(rows)
                        st.subheader("Prediction Table")
                        st.dataframe(df, use_container_width=True)
                    else:
                        st.code(text)

                except Exception as e:
                    st.error(f"Submission failed: {e}")
# ==================================================
# ABOUT PAGE
# ==================================================
elif section == "About":

    st.title("About")

    st.write("This tool performs protein characterization and isoform comparison using ProtParam-based analysis and UniProt sequence data.")
    st.write("This tool also predicts SubCellular Location by with recognition of signal peptide")

# ==================================================
# TEAM PAGE
# ==================================================
elif section == "Team":

    st.markdown("""
    <style>
    .team-header {
        display: flex;
        align-items: center;
        gap: 25px;
        margin-top: 10px;
        margin-bottom: 35px;
    }

    .team-header img {
        width: 110px;
        height: auto;
    }

    .team-heading-text {
        text-align: left;
    }

    .team-heading-text h1 {
        margin: 0;
        color: white !important;
        font-size: 40px;
        font-weight: 700;
    }

    .team-heading-text h2 {
        margin: 6px 0 0 0;
        color: white !important;
        font-size: 34px;
        font-weight: 600;
    }

    .team-heading-text h3 {
        margin: 5px 0 0 0;
        color: #c9d1d9 !important;
        font-size: 18px;
        font-weight: 400;
    }

    .team-heading-text h4 {
        margin: 8px 0 0 0;
        color: white !important;
        font-size: 18px;
        font-weight: 500;
    }

    .team-container {
        display: flex;
        justify-content: center;
        gap: 40px;
        flex-wrap: wrap;
        margin-top: 20px;
    }

    .team-card {
        background-color: #161b22;
        padding: 22px 28px;
        border-radius: 12px;
        text-align: center;
        width: 260px;
        border: 1px solid #30363d;
        box-shadow: 0 4px 15px rgba(0,0,0,0.25);
    }

    .team-card h3 {
        color: #58a6ff !important;
        margin-bottom: 12px;
        font-size: 20px;
    }

    .team-card p,
    .team-card a {
        color: white !important;
        font-size: 14px;
        text-decoration: none;
        margin: 0;
    }

    .team-card a:hover {
        color: #58a6ff !important;
        text-decoration: underline;
    }
    </style>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="team-header">

    <img src="https://despu.edu.in/assets/front/img/logo/xdes-white-logo.png.pagespeed.ic.ufvzDTbluy.png">

    <div class="team-heading-text">
        <h2>DES Pune University</h2>
        <h3>Department of Life Science</h3>
        <h4>Team</h4>
    </div>

    </div>

    <div class="team-container">

    <div class="team-card">
        <h3>Aakash Roy</h3>
        <p>Email:<br>
        <a href="mailto:royaakash2002@gmail.com">
        royaakash2002@gmail.com
        </a></p>
    </div>

    <div class="team-card">
        <h3>Snehal Rajkumar Yadav</h3>
        <p>Email:<br>
        <a href="mailto:3522511025@despu.edu.in">
        3522511025@despu.edu.in
        </a></p>
    </div>
    
    <div class="team-card">
        <h3>Harshada Risbud</h3>
        <p>Email:<br>
        <a href="mailto:3522511016@despu.edu.in">
        3522511016@despu.edu.in
        </a></p>
    </div>
    </div>
    """, unsafe_allow_html=True)
