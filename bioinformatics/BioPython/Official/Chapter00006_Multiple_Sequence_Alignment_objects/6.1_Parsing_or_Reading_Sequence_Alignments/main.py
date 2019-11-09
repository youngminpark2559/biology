

# We have two functions for reading in sequence alignments, Bio.AlignIO.read() and Bio.AlignIO.parse() which following the convention introduced in Bio.SeqIO are for files containing one or multiple alignments respectively.

# reading in sequence alignments
#   Bio.AlignIO.read(files containing one alignment)
#   Bio.AlignIO.parse(files containing multiple alignment)

# ================================================================================
# Using Bio.AlignIO.parse() will return an iterator which gives MultipleSeqAlignment objects. 

# iterator which gives MultipleSeqAlignment objects = Bio.AlignIO.parse(files which contain only a single alignment)

# ================================================================================
# Iterators are typically used in a for loop. 

# ================================================================================
# Examples of situations where you will have multiple different alignments include resampled alignments from the PHYLIP tool seqboot, or multiple pairwise alignments from the EMBOSS tools water or needle, or Bill Pearson’s FASTA tools.

# multiple different alignments
#   resampled alignments from the PHYLIP tool seqboot
#   multiple pairwise alignments from the EMBOSS tools water or needle
#   Bill Pearson’s FASTA tools

# ================================================================================
# However, in many situations you will be dealing with files which contain only a single alignment. In this case, you should use the Bio.AlignIO.read() function which returns a single MultipleSeqAlignment object.

# files which contain only a single alignment: general case
# single MultipleSeqAlignment object = Bio.AlignIO.read(files which contain only a single alignment)

# ================================================================================
# Both functions expect two mandatory arguments:

# ================================================================================
# The first argument is a handle to read the data from, typically an open file (see Section 24.1), or a filename.

# ================================================================================
# The second argument is a lower case string specifying the alignment format. 

# ================================================================================
# As in Bio.SeqIO we don’t try and guess the file format for you! See http://biopython.org/wiki/AlignIO for a full listing of supported formats.

# ================================================================================
# There is also an optional seq_count argument which is discussed in Section 6.1.3 below for dealing with ambiguous file formats which may contain more than one alignment.

# ================================================================================
# A further optional alphabet argument allowing you to specify the expected alphabet. 

# ================================================================================
# This can be useful as many alignment file formats do not explicitly label the sequences as RNA, DNA or protein – which means Bio.AlignIO will default to using a generic alphabet.

