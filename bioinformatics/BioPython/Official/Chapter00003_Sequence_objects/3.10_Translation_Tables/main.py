# conda activate py36gputorch100 && \
# cd /home/young/TEMP_GITHUB/biology/bioinformatics/BioPython/Official/Chapter00003_Sequence_objects/3.10_Translation_Tables && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
from Bio.Data import CodonTable

# ================================================================================
# In the previous sections we talked about the Seq object translation method (and mentioned the equivalent function in the Bio.Seq module – see Section 3.14). 

# ================================================================================
# Internally these use codon table objects derived from the NCBI information at ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt, also shown on https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi in a much more readable layout.

# ================================================================================
# As before, let’s just focus on two choices: the Standard translation table, and the translation table for Vertebrate Mitochondrial DNA.

# ================================================================================
translation_table_for_general_DNA = CodonTable.unambiguous_dna_by_name["Standard"]
translation_table_for_Mitochondrial_DNA = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

# ================================================================================
# Alternatively, these tables are labeled with ID numbers 1 and 2, respectively:

# ================================================================================
translation_table_for_general_DNA = CodonTable.unambiguous_dna_by_id[1]
translation_table_for_Mitochondrial_DNA = CodonTable.unambiguous_dna_by_id[2]

# ================================================================================
# You can compare the actual tables visually by printing them:

# ================================================================================
# print(translation_table_for_general_DNA)

# Table 1 Standard, SGC0

#   |  T      |  C      |  A      |  G      |
# --+---------+---------+---------+---------+--
# T | TTT F   | TCT S   | TAT Y   | TGT C   | T
# T | TTC F   | TCC S   | TAC Y   | TGC C   | C
# T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
# T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
# --+---------+---------+---------+---------+--
# C | CTT L   | CCT P   | CAT H   | CGT R   | T
# C | CTC L   | CCC P   | CAC H   | CGC R   | C
# C | CTA L   | CCA P   | CAA Q   | CGA R   | A
# C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
# --+---------+---------+---------+---------+--
# A | ATT I   | ACT T   | AAT N   | AGT S   | T
# A | ATC I   | ACC T   | AAC N   | AGC S   | C
# A | ATA I   | ACA T   | AAA K   | AGA R   | A
# A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
# --+---------+---------+---------+---------+--
# G | GTT V   | GCT A   | GAT D   | GGT G   | T
# G | GTC V   | GCC A   | GAC D   | GGC G   | C
# G | GTA V   | GCA A   | GAA E   | GGA G   | A
# G | GTG V   | GCG A   | GAG E   | GGG G   | G
# --+---------+---------+---------+---------+--

# ================================================================================
# print(translation_table_for_Mitochondrial_DNA)
# Table 2 Vertebrate Mitochondrial, SGC1

#   |  T      |  C      |  A      |  G      |
# --+---------+---------+---------+---------+--
# T | TTT F   | TCT S   | TAT Y   | TGT C   | T
# T | TTC F   | TCC S   | TAC Y   | TGC C   | C
# T | TTA L   | TCA S   | TAA Stop| TGA W   | A
# T | TTG L   | TCG S   | TAG Stop| TGG W   | G
# --+---------+---------+---------+---------+--
# C | CTT L   | CCT P   | CAT H   | CGT R   | T
# C | CTC L   | CCC P   | CAC H   | CGC R   | C
# C | CTA L   | CCA P   | CAA Q   | CGA R   | A
# C | CTG L   | CCG P   | CAG Q   | CGG R   | G
# --+---------+---------+---------+---------+--
# A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
# A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
# A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
# A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
# --+---------+---------+---------+---------+--
# G | GTT V   | GCT A   | GAT D   | GGT G   | T
# G | GTC V   | GCC A   | GAC D   | GGC G   | C
# G | GTA V   | GCA A   | GAA E   | GGA G   | A
# G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
# --+---------+---------+---------+---------+--

# ================================================================================
# You may find these following properties useful – for example if you are trying to do your own gene finding:

stop_codons_from_translation_table_for_Mitochondrial_DNA=translation_table_for_Mitochondrial_DNA.stop_codons
print("stop_codons_from_translation_table_for_Mitochondrial_DNA",stop_codons_from_translation_table_for_Mitochondrial_DNA)
# ['TAA', 'TAG', 'AGA', 'AGG']

start_codons_from_translation_table_for_Mitochondrial_DNA=translation_table_for_Mitochondrial_DNA.start_codons
print("start_codons_from_translation_table_for_Mitochondrial_DNA",start_codons_from_translation_table_for_Mitochondrial_DNA)
# ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']

translation_table_for_Mitochondrial_DNA.forward_table["ACG"]
# 'T'