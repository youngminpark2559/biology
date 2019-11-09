
# ================================================================================
# @@ NCBI BLAST+ (written in C++) was first released in 2009 as a replacement for the original NCBI “legacy” BLAST (written in C) which is no longer being updated.

# ================================================================================
# There were a lot of changes – the old version had a single core command line tool blastall which covered multiple different BLAST search types (which are now separate commands in BLAST+), and all the command line options were renamed.

# ================================================================================
# @@ Biopython’s wrappers for the NCBI “legacy” BLAST tools have been deprecated and will be removed in a future release.

# ================================================================================
# @@ To try to avoid confusion, we do not cover calling these old tools from Biopython in this tutorial.

# ================================================================================
# @@ You may also come across Washington University BLAST (WU-BLAST), and its successor, Advanced Biocomputing BLAST (AB-BLAST, released in 2009, not free/open source).

# ================================================================================
# @@ These packages include the command line tools wu-blastall and ab-blastall, which mimicked blastall from the NCBI “legacy” BLAST suite.

# ================================================================================
# @@ Biopython does not currently provide wrappers for calling these tools, but should be able to parse any NCBI compatible output from them.