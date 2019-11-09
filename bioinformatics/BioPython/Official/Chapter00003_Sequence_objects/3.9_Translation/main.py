# conda activate py36gputorch100 && \
# cd /home/young/TEMP_GITHUB/biology/bioinformatics/BioPython/Official/Chapter00003_Sequence_objects/3.9_Translation && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

# ================================================================================
# Sticking with the same example discussed in the transcription section above, now let’s translate this mRNA into the corresponding protein sequence - again taking advantage of one of the Seq object’s biological methods:

# mRNA -> translate -> corresponding protein sequence

# ================================================================================
# c messenger_rna: create mRNA sequence
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", IUPAC.unambiguous_rna)
# print("messenger_rna",messenger_rna)
# AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

protein_seq_which_is_translated_from_mRNA = messenger_rna.translate()
# print("protein_seq_which_is_translated_from_mRNA",protein_seq_which_is_translated_from_mRNA)
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

# ================================================================================
# You can also translate directly from the coding strand DNA sequence:

# coding strand DNA sequence -> directly translate -> corresponding protein sequence

# ================================================================================
# c coding_dna: DNA sequence (coding strand)
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
# print("coding_dna",coding_dna)
# Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

protein_seq_which_is_translated_from_DNA_sequence_coding_strand = coding_dna.translate()
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand",protein_seq_which_is_translated_from_DNA_sequence_coding_strand)
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

# ================================================================================
# You should notice in the above protein sequences that in addition to the end stop character, there is an internal stop as well. 

# ================================================================================
# This was a deliberate choice of example, as it gives an excuse to talk about some optional arguments, including different translation tables (Genetic Codes).

# ================================================================================
# The translation tables available in Biopython are based on those from the NCBI (see the next section of this tutorial). 
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

# ================================================================================
# By default, translation will use the standard genetic code (NCBI table id 1). 

# ================================================================================
# Suppose we are dealing with a mitochondrial sequence. 

# ================================================================================
# We need to tell the translation function to use the relevant genetic code instead:

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_Vertebrate_Mitochondrial_table=coding_dna.translate(table="Vertebrate Mitochondrial")
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_Vertebrate_Mitochondrial_table",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_Vertebrate_Mitochondrial_table)
# Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

# ================================================================================
# You can also specify the table using the NCBI table number which is shorter, and often included in the feature annotation of GenBank files:

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2=coding_dna.translate(table=2)
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2)
# Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

# ================================================================================
# Now, you may want to translate the nucleotides up to the first in frame stop codon, and then stop (as happens in nature):

protein_seq_which_is_translated_from_DNA_sequence_coding_strand = coding_dna.translate()
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand",protein_seq_which_is_translated_from_DNA_sequence_coding_strand)
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_translation_stops_at_stop_codon = coding_dna.translate(to_stop=True)
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_translation_stops_at_stop_codon",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_translation_stops_at_stop_codon)
# Seq('MAIVMGR', IUPACProtein())

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2 = coding_dna.translate(table=2)
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2)
# MAIVMGRWKGAR*

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_translation_stops_at_stop_codon = coding_dna.translate(table=2, to_stop=True)
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_translation_stops_at_stop_codon",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_translation_stops_at_stop_codon)
# MAIVMGRWKGAR

# ================================================================================
# Notice that when you use the to_stop argument, the stop codon itself is not translated - and the stop symbol is not included at the end of your protein sequence.

# ================================================================================
# You can even specify the stop symbol if you don’t like the default asterisk:

protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_using_custom_stop_codon_symbol = coding_dna.translate(table=2, stop_symbol="@")
# print("protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_using_custom_stop_codon_symbol",protein_seq_which_is_translated_from_DNA_sequence_coding_strand_using_table_2_using_custom_stop_codon_symbol)
# Seq('MAIVMGRWKGAR@', HasStopCodon(IUPACProtein(), '@'))

# ================================================================================
# Now, suppose you have a complete coding sequence CDS, which is to say a nucleotide sequence (e.g. mRNA – after any splicing) which is a whole number of codons (i.e. the length is a multiple of three), commences with a start codon, ends with a stop codon, and has no internal in-frame stop codons. 

# complete coding sequence: 
#   nucleotide sequence
#   mRNA – after any splicing
#   whole number of codons
#   commences with a start codon, ends with a stop codon
#   has no internal in-frame stop codons

# ================================================================================
# In general, given a complete CDS, the default translate method will do what you want (perhaps with the to_stop option). 

# Works well: default translate(default translate)

# ================================================================================
# However, what if your sequence uses a non-standard start codon? 

# sequence: 
# non-standard start codon

# ================================================================================
# This happens a lot in bacteria – for example the gene yaaX in E. coli K12:

# ================================================================================
# c gene: complete coding sequence
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
           "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
           "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
           "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
           "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
           generic_dna)


protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table = gene.translate(table="Bacterial")
# print("protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table",protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table)
# Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*', HasStopCodon(ExtendedIUPACProtein(), '*')

protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table_stops_at_stop_codon = gene.translate(table="Bacterial", to_stop=True)
# print("protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table_stops_at_stop_codon",protein_which_is_translated_from_complete_coding_sequence_using_Bacterial_table_stops_at_stop_codon)
# Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR', ExtendedIUPACProtein())

# ================================================================================
# In the bacterial genetic code, GTG is a valid start codon, and while it does normally encode Valine, if used as a start codon it should be translated as methionine. This happens if you tell Biopython your sequence is a complete CDS:

# bacterial genetic code
#   GTG is a valid start codon
#   GTG normally encodes Valine
#   if used as a start codon it should be translated as methionine

# ================================================================================
# protein_seq_which_is_translated_from_sequence_cds = gene.translate(table="Bacterial", cds=True)
# Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR', ExtendedIUPACProtein())

# ================================================================================
# In addition to telling Biopython to translate an alternative start codon as methionine, using this option also makes sure your sequence really is a valid CDS (you’ll get an exception if not).

# The example in Section 20.1.3 combines the Seq object’s translate method with Bio.SeqIO for sequence input/output.
