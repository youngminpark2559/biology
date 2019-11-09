# This chapter is about Multiple Sequence Alignments, by which we mean a collection of multiple sequences which have been aligned together – usually with the insertion of gap characters, and addition of leading or trailing gaps – such that all the sequence strings are the same length. 

# Multiple Sequence Alignments
#   collection of multiple sequences which have been aligned together
#   usually with the insertion of gap characters
#   leading gap or trailing gap
#   all the sequence strings are the same length. 

# ================================================================================
# Such an alignment can be regarded as a matrix of letters, where each row is held as a SeqRecord object internally.

# alignment
#   matrix of letters
#   each row is a SeqRecord object

# ================================================================================
# We will introduce the MultipleSeqAlignment object which holds this kind of data, and the Bio.AlignIO module for reading and writing them as various file formats (following the design of the Bio.SeqIO module from the previous chapter). 

# MultipleSeqAlignment object contains
#   matrix of letters
#   each row is a SeqRecord object

# Bio.AlignIO
#   reading and writing MultipleSeqAlignment objects as various file formats
  
# ================================================================================
# Note that both Bio.SeqIO and Bio.AlignIO can read and write sequence alignment files. 

# Bio.SeqIO / Bio.AlignIO
#   read and write sequence alignment files

# ================================================================================
# The appropriate choice will depend largely on what you want to do with the data.

# ================================================================================
# The final part of this chapter is about our command line wrappers for common multiple sequence alignment tools like ClustalW and MUSCLE.

