# conda activate py36gputorch100 && \
# cd /home/young/TEMP_GITHUB/biology/bioinformatics/BioPython/Official/Chapter00003_Sequence_objects/3.7_Nucleotide_sequences_and_reverse_complements && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# ================================================================================
# For nucleotide sequences, you can easily obtain the complement or reverse complement of a Seq object using its built-in methods:

# nucleotide sequences -> complement
# nucleotide sequences -> reverse complement

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
# print("my_seq",my_seq)
# GATCGATGGGCCTATATAGGATCGAAAATCGC
# print("my_seq",type(my_seq))
# <class 'Bio.Seq.Seq'>

# ================================================================================
complement_of_my_seq = my_seq.complement()
# print("complement_of_my_seq",complement_of_my_seq)
# CTAGCTACCCGGATATATCCTAGCTTTTAGCG
# print("complement_of_my_seq",type(complement_of_my_seq))
# <class 'Bio.Seq.Seq'>

# ================================================================================
reverse_complement_of_my_seq = my_seq.reverse_complement()
# print("reverse_complement_of_my_seq",reverse_complement_of_my_seq)
# GCGATTTTCGATCCTATATAGGCCCATCGATC
# print("reverse_complement_of_my_seq",type(reverse_complement_of_my_seq))
# <class 'Bio.Seq.Seq'>

# ================================================================================
# As mentioned earlier, an easy way to just reverse a Seq object (or a Python string) is slice it with -1 step:

reverse_of_my_seq = my_seq[::-1]
# print("reverse_of_my_seq",reverse_of_my_seq)
# CGCTAAAAGCTAGGATATATCCGGGTAGCTAG
# print("reverse_of_my_seq",type(reverse_of_my_seq))
# <class 'Bio.Seq.Seq'>

# ================================================================================
# In all of these operations, the alphabet property is maintained. 

# ================================================================================
# This is very useful in case you accidentally end up trying to do something weird like take the (reverse)complement of a protein sequence:

# ================================================================================
# protein_seq = Seq("EVRNAK", IUPAC.protein)
# complement_of_protein_seq = protein_seq.complement()
# print("complement_of_protein_seq",complement_of_protein_seq)
# Traceback (most recent call last):
#   File "main.py", line 53, in <module>
#     complement_of_protein_seq = protein_seq.complement()
#   File "/home/young/anaconda3/envs/py36gputorch100/lib/python3.6/site-packages/Bio/Seq.py", line 920, in complement
#     raise ValueError("Proteins do not have complements!")
# ValueError: Proteins do not have complements!

# ================================================================================
# The example in Section 5.5.3 combines the Seq object’s reverse complement method with Bio.SeqIO for sequence input/output.
