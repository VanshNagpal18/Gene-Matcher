import streamlit as st
import pandas as pd
import time

# -----------------------------------
# Gene Matching App using KMP & Boyer-Moore
# -----------------------------------

st.set_page_config(page_title="🧬 Gene Matcher", layout="wide")

st.title("🧬 Gene Matching Platform")
st.markdown("""
Compare *KMP* and *Boyer–Moore* algorithms for DNA sequence pattern matching.
Upload your DNA data or paste directly, choose algorithm, and analyze efficiently.
""")

# ---- Sidebar ----
st.sidebar.header("🔧 Settings")
algorithm_choice = st.sidebar.selectbox("Select Algorithm", ["KMP", "Boyer–Moore", "Both"])
st.sidebar.info("Developed by Vansh Nagpal | DAA Project")

# ---- DNA Sequence Input ----
uploaded_file = st.file_uploader("📤 Upload DNA Sequence File (TXT/FASTA)", type=["txt", "fasta"])

if uploaded_file:
    dna_sequence = uploaded_file.read().decode("utf-8").strip().replace("\n", "")
else:
    dna_sequence = st.text_area("Or Paste DNA Sequence Here:", height=150)

pattern = st.text_input("🔍 Enter Gene Pattern to Search:")

# ---- Helper Functions ----
def compute_lps(pattern):
    lps = [0] * len(pattern)
    length, i = 0, 1
    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(text, pattern):
    lps = compute_lps(pattern)
    i = j = 0
    matches = []
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == len(pattern):
            matches.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            j = lps[j - 1] if j != 0 else 0
            if j == 0: i += 1
    return matches

def bad_char_heuristic(pattern):
    bad_char = [-1] * 256
    for i in range(len(pattern)):
        bad_char[ord(pattern[i])] = i
    return bad_char

def boyer_moore_search(text, pattern):
    bad_char = bad_char_heuristic(pattern)
    m, n = len(pattern), len(text)
    matches = []
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            matches.append(s)
            s += (m - bad_char[ord(text[s + m])] if s + m < n else 1)
        else:
            s += max(1, j - bad_char[ord(text[s + j])])
    return matches

# ---- Main Logic ----
if st.button("🔬 Run Matching"):
    if dna_sequence and pattern:
        st.subheader("Results")

        results = []  # ✅ Proper indentation here

        if algorithm_choice in ["KMP", "Both"]:
            start = time.time()
            matches = kmp_search(dna_sequence, pattern)
            end = time.time()
            results.append(["KMP", len(matches), end - start])
            st.success(f"KMP found {len(matches)} matches in {end - start:.4f} seconds")

        if algorithm_choice in ["Boyer–Moore", "Both"]:
            start = time.time()
            matches = boyer_moore_search(dna_sequence, pattern)
            end = time.time()
            results.append(["Boyer–Moore", len(matches), end - start])
            st.success(f"Boyer–Moore found {len(matches)} matches in {end - start:.4f} seconds")

        # Display comparison table
        df = pd.DataFrame(results, columns=["Algorithm", "Matches Found", "Execution Time (s)"])
        st.table(df)

        # Save results
        csv_data = df.to_csv(index=False).encode("utf-8")
        st.download_button("💾 Download Results as CSV", csv_data, "gene_results.csv", "text/csv")

        # Generate shareable link
        st.markdown("🔗 *Shareable Results Link (copy manually):*")
        st.code(st.get_option("server.address") + f"?pattern={pattern}", language="text")

    else:
        st.warning("Please upload or paste a DNA sequence and enter a pattern.")
        
# ---- Footer ----
st.markdown("""
---
👩‍💻 *Project by Vansh Nagpal*  
Built using [Streamlit](https://streamlit.io/) for educational purposes.  
""")

