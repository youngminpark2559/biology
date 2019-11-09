# conda activate py36gputorch100 && \
# cd /mnt/external_disk/Code_projects/Bioinformatics/BioPython_official/Chapter00003_Sequence_objects/3.8_Transcription && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# ================================================================================
# @@ Before talking about transcription, I want to try to clarify the strand issue. 

# @@ Consider the following (made up) stretch of double stranded DNA which encodes a short peptide:

# ================================================================================
# DNA coding strand (aka Crick strand, strand +1)
# 5’	ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG	3’
#  	  |||||||||||||||||||||||||||||||||||||||
# 3’	TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC	5’
# DNA template strand (aka Watson strand, strand −1)
 
# |
# Transcription
# ↓

# 5’	AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG	3’
# Single stranded messenger RNA
 
# ================================================================================
# The actual biological transcription process works from the template strand, doing a reverse complement (TCAG → CUGA) to give the mRNA. 

# biological transcription process
# template strand
# reverse complement
# mRNA

# ================================================================================
# However, in Biopython and bioinformatics in general, we typically work directly with the coding strand because this means we can get the mRNA sequence just by switching T → U.

# ================================================================================
# Now let’s actually get down to doing a transcription in Biopython. 

# ================================================================================
# First, let’s create Seq objects for the coding and template DNA strands:

# ================================================================================
# c coding_dna: dna sequence
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
# print("coding_dna",coding_dna)
# ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

# c template_dna: template_dna which is generated from coding_dna via the process of reverse_complement
template_dna = coding_dna.reverse_complement()
# print("template_dna",template_dna)
# CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT

# ================================================================================
# These should match the figure above - remember by convention nucleotide sequences are normally read from the 5’ to 3’ direction, while in the figure the template strand is shown reversed.

# nucleotide sequences should be read from the 5’ to 3’

# ================================================================================
# Now let’s transcribe the coding strand into the corresponding mRNA, using the Seq object’s built in transcribe method:

# print("coding_dna",coding_dna)
# ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

# c messenger_rna: mRNA which is generated from DNA sequence (5' -> 3` coding strand)
messenger_rna = coding_dna.transcribe()
# print("messenger_rna",messenger_rna)
# AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

# ================================================================================
# As you can see, all this does is switch T → U, and adjust the alphabet.

# ================================================================================
# If you do want to do a true biological transcription starting with the template strand, then this becomes a two-step process:

# true biological transcription

# ================================================================================
# Step1: 
# template_dna -> reverse_complement() -> reverse_complement_of_template_dna
# Step2: 
# reverse_complement_of_template_dna -> transcribe()

mrna=template_dna.reverse_complement().transcribe()
# print("mrna",mrna)
# AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

# ================================================================================
# The Seq object also includes a back-transcription method for going from the mRNA to the coding strand of the DNA. 

# ================================================================================
# Again, this is a simple U → T substitution and associated change of alphabet:

# ================================================================================
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", IUPAC.unambiguous_rna)
# print("messenger_rna",messenger_rna)
# Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

# coding_sequence=messenger_rna.back_transcribe()
# print("coding_sequence",coding_sequence)
# Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

# ================================================================================
# Note: The Seq object’s transcribe and back_transcribe methods were added in Biopython 1.49. 

# ================================================================================
# For older releases you would have to use the Bio.Seq module’s functions instead, see Section 3.14.
