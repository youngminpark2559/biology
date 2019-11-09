# conda activate py36_network_medicine && \
# cd /mnt/external_disk/Code_projects/Bioinformatics/BioPython_official/Chapter00007_BLAST && \
# rm e.l && python 7.1_Running_BLAST_over_the_Internet.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
import time,timeit,datetime
from Bio.Blast import NCBIWWW
from Bio import SeqIO

# ================================================================================
# # \. 
# # .\n\n ================================================================================\n

# # ================================================================================
# 7.1  Running BLAST over the Internet

# # ================================================================================
# We use the function qblast() in the Bio.Blast.NCBIWWW module to call the online version of BLAST.

# call the online version of BLAST:
# Bio.Blast.NCBIWWW.qblast()

# # ================================================================================
# This has three non-optional arguments:

# # ================================================================================
# The first argument is the blast program to use for the search, as a lower case string.

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    )

# # ================================================================================
# The options and descriptions of the programs are available at https://blast.ncbi.nlm.nih.gov/Blast.cgi.

# # ================================================================================
# Currently qblast only works with blastn, blastp, blastx, tblast and tblastx.

# # ================================================================================
# The second argument specifies the databases to search against.

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    arg2=databases to search against
#    )

# # ================================================================================
# Again, the options for this are available on the NCBI Guide to BLAST ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf.

# # ================================================================================
# The third argument is a string containing your query sequence.

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    arg2=databases to search against,
#    arg3=string containing your query sequence)

# # ================================================================================
# This can either be the sequence itself, the sequence in fasta format, or an identifier like a GI number.

# arg3: "sequence itself", "sequence in fasta format", "identifier like a GI number"

# # ================================================================================
# The qblast function also take a number of other option arguments which are basically analogous to the different parameters you can set on the BLAST web page.

# # ================================================================================
# We’ll just highlight a few of them here:

# The argument url_base sets the base URL for running BLAST over the internet.

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    arg2=databases to search against,
#    arg3=string containing your query sequence,
#    url_base=base URL for running BLAST over the internet)

# # ================================================================================
# By default it connects to the NCBI, but one can use this to connect to an instance of NCBI BLAST running in the cloud.

# url_base: "NCBI", "instance of NCBI BLAST running in the cloud"

# # ================================================================================
# @@ Please refer to the documentation for the qblast function for further details.

# # ================================================================================
# The qblast function can return the BLAST results in various formats, which you can choose with the optional format_type keyword: "HTML", "Text", "ASN.1", or "XML".

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    arg2=databases to search against,
#    arg3=string containing your query sequence,
#    url_base=base URL for running BLAST over the internet,
#    format_type=BLAST results in various formats)

# format_type: "HTML", "Text", "ASN.1", "XML"

# # ================================================================================
# @@ The default is "XML", as that is the format expected by the parser, described in section 7.3 below.

# # ================================================================================
# The argument expect sets the expectation or e-value threshold.

# Bio.Blast.NCBIWWW.qblast(
#    arg1=blast program to use for the search,
#    arg2=databases to search against,
#    arg3=string containing your query sequence,
#    url_base=base URL for running BLAST over the internet,
#    format_type=BLAST results in various formats,
#    expectation=e-value threshold)

# # ================================================================================
# @@ For more about the optional BLAST arguments, we refer you to the NCBI’s own documentation, or that built into Biopython:

# # ================================================================================
# from Bio.Blast import NCBIWWW
# help(NCBIWWW.qblast)
# ...

# ================================================================================
# Note that the default settings on the NCBI BLAST website are not quite the same as the defaults on QBLAST. 

# ================================================================================
# If you get different results, you’ll need to check the parameters (e.g., the expectation value threshold and the gap values).

# ================================================================================
# For example, if you have a nucleotide sequence you want to search against the nucleotide database (nt) using BLASTN, and you know the GI number of your query sequence, you can use:

# "nucleotide sequence" you want to search against the "nucleotide database (nt)" using "BLASTN"
# "GI number" of your query sequence

# ================================================================================
# start=timeit.default_timer()
# result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
# print("result_handle",result_handle)
# <_io.StringIO object at 0x7fcd8b348318>

# stop=timeit.default_timer()
# took_time_sec=stop-start
# took_time_min=str(datetime.timedelta(seconds=took_time_sec))
# print('took_time_min',took_time_min)
# took_time_min 0:01:13.794526

# ================================================================================
# Alternatively, if we have our "query sequence" already in a "FASTA formatted file", we just need to "open the file" and "read in this record as a string", and use that as the "query argument":

# m_cold.fasta: FASTA formatted file in your local PC
# fasta_string = open("m_cold.fasta").read()

# blastn: the kind of BLAST program
# nt: the kind of database (nucleotide database)
# result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

# ================================================================================
# We could also have read in the FASTA file as a SeqRecord and then supplied just the sequence itself:

# Pass m_cold.fasta file to SeqIO.read() with specifying that it's fasta formatted file
# record = SeqIO.read("m_cold.fasta", format="fasta")

# record.seq: Extract sequence from record
# result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# ================================================================================
# Supplying just the sequence means that BLAST will assign an identifier for your sequence automatically. 

# You might prefer to use the SeqRecord object’s format method to make a FASTA string (which will include the existing identifier):

# record = SeqIO.read("m_cold.fasta", format="fasta")
# result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

# ================================================================================
# This approach makes more sense if you have your sequence(s) in a non-FASTA file format which you can extract using Bio.SeqIO (see Chapter 5).

# ================================================================================
# Whatever arguments you give the qblast() function, you should get back your results in a handle object (by default in XML format).

# ================================================================================
# The next step would be to parse the XML output into Python objects representing the search results (Section 7.3), but you might want to save a local copy of the output file first.

# saved file = save(local copy of the output file)

# Python objects = parse(the XML output)

# ================================================================================
# I find this especially useful when debugging my code that extracts info from the BLAST results (because re-running the online search is slow and wastes the NCBI computer time) 

# ================================================================================
# We need to be a bit careful since we can use result_handle.read() to read the BLAST output only once – calling result_handle.read() again returns an empty string.

# ================================================================================
# Open my_blast.xml file in write mode with creating file stream named out_handle
with open("my_blast.xml", "w") as out_handle:
   read_result=result_handle.read()
   out_handle.write(read_result)

# save read_result into my_blast.xml
result_handle.close()

# ================================================================================
# After doing this, the results are in the file my_blast.xml and the original handle has had all its data extracted (so we closed it). 

# ================================================================================
# However, the parse function of the BLAST parser (described in 7.3) takes a file-handle-like object, so we can just open the saved file for input:
# result_handle = open("my_blast.xml")

# ================================================================================
# Now that we’ve got the BLAST results back into a handle again, we are ready to do something with them, so this leads us right into the parsing section (see Section 7.3 below). 

# ================================================================================
# You may want to jump ahead to that now ….
