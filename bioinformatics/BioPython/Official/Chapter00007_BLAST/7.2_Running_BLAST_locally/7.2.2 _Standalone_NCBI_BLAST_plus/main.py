
# ================================================================================
# @@ The “new” NCBI BLAST+ suite was released in 2009. 

# ================================================================================
# @@ This replaces the old NCBI “legacy” BLAST package (see below).

# ================================================================================
# @@ This section will show briefly how to use these tools from within Python. 

# ================================================================================
# @@ If you have already read or tried the alignment tool examples in Section 6.4 this should all seem quite straightforward. 

# ================================================================================
# First, we construct a command line string (as you would type in at the command line prompt if running standalone BLAST by hand). 

# ================================================================================
# Then we can execute this command from within Python.

# For example, taking a FASTA file of gene nucleotide sequences, you might want to run a BLASTX (translation) search against the non-redundant (NR) protein database. 

# ================================================================================
# Assuming you (or your systems administrator) has downloaded and installed the NR database, you might run:

# FASTA file of gene nucleotide sequences
# run a BLASTX (translation) search
# against the non-redundant (NR) protein database
# NR database in your local PC

# ================================================================================
# expectation cut-off value of 0.001
# $ blastx -query opuntia.fasta -db nr -out opuntia.xml -evalue 0.001 -outfmt 5

# ================================================================================
# This should run BLASTX against the NR database, using an expectation cut-off value of 0.001 and produce XML output to the specified file (which we can then parse). 

# On my computer this takes about six minutes - a good reason to save the output to a file so you can repeat any analysis as needed.

# ================================================================================
# From within Biopython we can use the NCBI BLASTX wrapper from the Bio.Blast.Applications module to build the command line string, and run it:

# from Bio.Blast.Applications import NcbiblastxCommandline
# help(NcbiblastxCommandline)

# ================================================================================
# blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")
# blastx_cline
# NcbiblastxCommandline(cmd='blastx', out='opuntia.xml', outfmt=5, query='opuntia.fasta', db='nr', evalue=0.001)

# print(blastx_cline)
# blastx -out opuntia.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001

# stdout, stderr = blastx_cline()

# ================================================================================
# In this example there shouldn’t be any output from BLASTX to the terminal, so stdout and stderr should be empty. 

# ================================================================================
# You may want to check the output file opuntia.xml has been created.

# ================================================================================
# As you may recall from earlier examples in the tutorial, the opuntia.fasta contains seven sequences, so the BLAST XML output should contain multiple results. 

# ================================================================================
# Therefore use Bio.Blast.NCBIXML.parse() to parse it as described below in Section 7.3.

