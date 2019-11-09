# conda activate py36gputorch100 && \
# cd /mnt/external_disk/Code_projects/Bioinformatics/BioPython_official/Chapter00007_BLAST/7.3_Parsing_BLAST_output && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# ================================================================================
# As mentioned above, BLAST can generate output in various formats, such as XML, HTML, and plain text.

# result = BLAST(sequence)
# result: "XML", "HTML", "plain text"

# ================================================================================
# @@ Originally, Biopython had parsers for BLAST plain text and HTML output, as these were the only output formats offered at the time.

# ================================================================================
# @@ Unfortunately, the BLAST output in these formats kept changing, each time breaking the Biopython parsers.

# ================================================================================
# @@ Our HTML BLAST parser has been removed, while the deprectaed plain text BLAST parser is now only available via Bio.SearchIO.

# ================================================================================
# @@ Use it at your own risk, it may or may not work, depending on which BLAST version you’re using.

# ================================================================================
# @@ As keeping up with changes in BLAST became a hopeless endeavor, especially with users running different BLAST versions, we now recommend to parse the output in XML format, which can be generated by recent versions of BLAST.

# ================================================================================
# @@ Not only is the XML output more stable than the plain text and HTML output, it is also much easier to parse automatically, making Biopython a whole lot more stable.

# ================================================================================
# @@ You can get BLAST output in XML format in various ways.

# ================================================================================
# @@ For the parser, it doesn’t matter how the output was generated, as long as it is in the XML format.

# ================================================================================
# @@ You can use Biopython to run BLAST over the internet, as described in section 7.1.

# ================================================================================
# @@ You can use Biopython to run BLAST locally, as described in section 7.2.

# ================================================================================
# @@ You can do the BLAST search yourself on the NCBI site through your web browser, and then save the results.

# ================================================================================
# @@ You need to choose XML as the format in which to receive the results, and save the final BLAST page you get (you know, the one with all of the interesting results!) to a file.

# ================================================================================
# @@ You can also run BLAST locally without using Biopython, and save the output in a file.

# ================================================================================
# @@ Again, you need to choose XML as the format in which to receive the results.

# ================================================================================
# @@ The important point is that you do not have to use Biopython scripts to fetch the data in order to be able to parse it.

# ================================================================================
# @@ Doing things in one of these ways, you then need to get a handle to the results.

# ================================================================================
# @@ In Python, a handle is just a nice general way of describing input to any info source so that the info can be retrieved using read() and readline() functions (see Section 24.1).

# ================================================================================
# @@ If you followed the code above for interacting with BLAST through a script, then you already have result_handle, the handle to the BLAST results.

# ================================================================================
# For example, using a GI number to do an online search:

# blastn: the kind of BLAST program
# nt: database
# 8332116: GI number
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
# print("result_handle",result_handle)
# <_io.StringIO object at 0x7fbbf4f53558>

# print(result_handle.read())
# <?xml version="1.0"?>
# <!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
# <BlastOutput>
#   <BlastOutput_program>blastn</BlastOutput_program>
#   <BlastOutput_version>BLASTN 2.10.0+</BlastOutput_version>
#   <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
#   <BlastOutput_db>nt_v5</BlastOutput_db>
#   <BlastOutput_query-ID>BE037100.1</BlastOutput_query-ID>
#   <BlastOutput_query-def>MP14H09 MP Mesembryanthemum crystallinum cDNA 5&apos; similar to cold acclimation protein, mRNA sequence</BlastOutput_query-def>
#   <BlastOutput_query-len>1111</BlastOutput_query-len>
#   <BlastOutput_param>
#     <Parameters>
#       <Parameters_expect>10</Parameters_expect>
#       <Parameters_sc-match>2</Parameters_sc-match>
#       <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
#       <Parameters_gap-open>5</Parameters_gap-open>
#       <Parameters_gap-extend>2</Parameters_gap-extend>
#       <Parameters_filter>L;m;</Parameters_filter>
#     </Parameters>
#   </BlastOutput_param>
# <BlastOutput_iterations>
# <Iteration>
#   <Iteration_iter-num>1</Iteration_iter-num>
#   <Iteration_query-ID>BE037100.1</Iteration_query-ID>
#   <Iteration_query-def>MP14H09 MP Mesembryanthemum crystallinum cDNA 5&apos; similar to cold acclimation protein, mRNA sequence</Iteration_query-def>
#   <Iteration_query-len>1111</Iteration_query-len>
# <Iteration_hits>
# <Hit>
#   <Hit_num>1</Hit_num>
#   <Hit_id>gi|1219041180|ref|XM_021875076.1|</Hit_id>
#   <Hit_def>PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA</Hit_def>
#   <Hit_accession>XM_021875076</Hit_accession>
#   <Hit_len>1173</Hit_len>
#   <Hit_hsps>
#     <Hsp>
#       <Hsp_num>1</Hsp_num>
#       <Hsp_bit-score>435.898</Hsp_bit-score>
#       <Hsp_score>482</Hsp_score>
#       <Hsp_evalue>1.65665e-117</Hsp_evalue>
#       <Hsp_query-from>59</Hsp_query-from>
#       <Hsp_query-to>678</Hsp_query-to>
#       <Hsp_hit-from>278</Hsp_hit-from>
#       <Hsp_hit-to>901</Hsp_hit-to>
#       <Hsp_query-frame>1</Hsp_query-frame>
#       <Hsp_hit-frame>1</Hsp_hit-frame>
#       <Hsp_identity>473</Hsp_identity>
#       <Hsp_positive>473</Hsp_positive>
#       <Hsp_gaps>4</Hsp_gaps>
#       <Hsp_align-len>624</Hsp_align-len>
#       <Hsp_qseq>ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGT--CATTACGGG-TTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCACTTTTAGTCCCTTACATATTCCTCAGTCTTCCATCCGGGCCATTTCATCTGTTCAGAGGCGAGGTCGGGAAATGGATTGCCATCATTGCAGTCGTGTTAAGGCTGTTCTTCAACCGGCATTTCCCAGTTTGGCTGGAAATGCCTGGATCGTTGATACTCCTCCTGGTGGTGGCACCAGACTTCTTTACACACAAAGTGAAGGAGAGCTGGATCGGAATTGCAATTATGATAGCGATAGGGTGTCACCTGATGCAAGAACATATCAGAGCCACTGGTGGCTTTTGGAATTCCTTCACACAGAGCCACGGAACTTTTAACACAATTGGGCTTATCCTTCTACTGGCTTACCCTGTCT-GTTTATGGTCATCTTCATGATGTA</Hsp_qseq>
#       <Hsp_hseq>ACCGAAAATGGGCAGAGGAGTGAATTATATGGCAATGACACCTGAGCAACTAGCCGCGGCCAATTTGATCAACTCCGACATCAATGAGCTCAAGATCGTTGTGATGACACTCATTCATGATGCTTCTAGACTCGGCGGCACCTCAGGATTTGGAACTCATTTTCTTAGATGGCTAGCCTCTCTTGCTGCTATTTACTTGTTGATCCTGGATCGCACAAATTGGAGAACCAACATGCTCACATCACTCTTAGTACCATACATATTCCTCAGTCTTCCTTCTGGCCCTTTTTACCTTCTTAGGGGTGAGGTTGGGAAATGGATTGCTTTTGTCGCGGTTGTGCTAAGGCTATTCTTCCACCGCCGCTTCCCAGAATGGTTAGAGATGCCAGGATCACTGATACTATTGTTGGTGGTAGCTCCAGAATTGCTAGCACACAAATTAAAGGATAGTTGGATGGGAGTTGTAATTCTGTTAATCATAGGGTGTTATTTGCTGCAAGAACATATCAGGGCAACTGGTGGTTTAAGAAATTCGTTTACTCAAAGCCATGGAATTTCCTATACGATTGGGCTGCTTCTCTTATTGGCTTACCCAATTTGGTCCATGGTTATTTTCATGATTTA</Hsp_hseq>
#       <Hsp_midline>|| ||||||||| |||| | |||| ||  |||| |||| | |||| ||| | |||| ||| ||| ||||| | ||||| ||||||||||| || || |     ||||  |||||  ||||||||  ||  |||||   ||   | || ||||| |||||||| || | |||||| ||||  ||||| |||||||||||||||||  ||||||| ||||| |||||||||||||||||||| ||||| ||||| || |||||||||||||||||||| || || || ||| | ||  | || || ||||| ||||||||||||||  |  | || || ||| ||||||| |||||| |||| |  |||||||  ||| | || ||||| |||||  |||||||  |  ||||||| || ||||| ||  |  |||||||| | ||||| || ||||| ||| ||| |||| || ||   ||||||||| |  || |||||||||||||||| || |||||||| ||  | ||||| || || || ||||| |||| ||   | || ||||||||  | ||  || ||||||||||  | | ||  ||||| || |||||||| ||</Hsp_midline>
#     </Hsp>
#   </Hit_hsps>
# </Hit>
# <Hit>
#   <Hit_num>2</Hit_num>
#   <Hit_id>gi|1226796956|ref|XM_021992092.1|</Hit_id>
#   <Hit_def>PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA</Hit_def>
#   <Hit_accession>XM_021992092</Hit_accession>
#   <Hit_len>672</Hit_len>
#   <Hit_hsps>
#     <Hsp>
#       <Hsp_num>1</Hsp_num>
#       <Hsp_bit-score>423.275</Hsp_bit-score>
#       <Hsp_score>468</Hsp_score>
#       <Hsp_evalue>1.04546e-113</Hsp_evalue>
#       <Hsp_query-from>63</Hsp_query-from>
#       <Hsp_query-to>649</Hsp_query-to>
#       <Hsp_hit-from>11</Hsp_hit-from>
#       <Hsp_hit-to>600</Hsp_hit-to>
#       <Hsp_query-frame>1</Hsp_query-frame>
#       <Hsp_hit-frame>1</Hsp_hit-frame>
#       <Hsp_identity>448</Hsp_identity>
#       <Hsp_positive>448</Hsp_positive>
#       <Hsp_gaps>3</Hsp_gaps>
#       <Hsp_align-len>590</Hsp_align-len>
#       <Hsp_qseq>AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGG---TCATTACGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCACTTTTAGTCCCTTACATATTCCTCAGTCTTCCATCCGGGCCATTTCATCTGTTCAGAGGCGAGGTCGGGAAATGGATTGCCATCATTGCAGTCGTGTTAAGGCTGTTCTTCAACCGGCATTTCCCAGTTTGGCTGGAAATGCCTGGATCGTTGATACTCCTCCTGGTGGTGGCACCAGACTTCTTTACACACAAAGTGAAGGAGAGCTGGATCGGAATTGCAATTATGATAGCGATAGGGTGTCACCTGATGCAAGAACATATCAGAGCCACTGGTGGCTTTTGGAATTCCTTCACACAGAGCCACGGAACTTTTAACACAATTGGGCTTATCCTTCTACTGGCTTACCC</Hsp_qseq>
#       <Hsp_hseq>AAAATGGGTAGACGAATGGATTATTTGGCGATGAAAACCGAGCAATTAGCCGCGGCCAATTTGATCGATTCCGATATCAACGAGCTGAAGATCGCCGTGATGGCGCTCGTTCATGATACTACTACGCTCGGCGGTCAATCGGGATTTGGAACTCATTTTCTCCAATGGCTCGCCTCATTTTCTGCTATTTACTTGTTAATCCTTGATCGAACACATTGGAGAAGCAACATGCTTACTTCACTTTTAGTACCATACATTTTCCTAAGTCTCCCATCTGGCCCCTTTCACCTTTTAAGAGGTGAGGTTGGGAAATGGATTGCTTTTGTCTCGGTTGTGCTAAGGTTATTCTTCCACCGCAGTTTCCCAGAATGGTTGGAAATGCCAGTATGTTTGATACTATTATTGGTGGTAGCTCCAGAAATGCTTGCAATATCAATGAAAGAGAGTTGGATGGGAGTTGTAGTTGTGTTAATCATAGGATGTTACCTTCTACAAGAGCATATTAGGGCAACTGGTGGTTTAAGGAATTCTTTCACACAAAGACATGGGATTTCCAACACAATTGGGCTTCTTCTCTTGTTGGCTTACCC</Hsp_hseq>
#       <Hsp_midline>|||||||| |||  |||| | || ||||| |||||||| || ||||| |||| ||| ||| ||||||||||||||||||| ||||| || || ||    |||  |||| |  ||||| ||| || ||||||   ||| |  || ||||| |||||||| ||| ||||||||||||   || | |||||||||||||| ||  | ||||||||| | ||||||| ||||||||| || ||||||||||| || ||||| ||||| ||||| ||||| || || ||||| || || ||||| ||||| ||||||||||||||  |  |  | || ||| ||||| | |||||| ||||   ||||||||  ||| |||||||||| | ||  ||||||||  |  ||||||| || |||||  |  || ||     | |||| ||||| ||||| ||| ||| | || || ||   ||||| ||| ||||  | ||||| ||||| || || |||||||| ||  ||||||| |||||||| || || || | ||  ||||||||||||||| | ||  |  ||||||||||</Hsp_midline>
#     </Hsp>
#   </Hit_hsps>
# ...

# ================================================================================
# If instead you ran BLAST some other way, and have the BLAST output (in XML format) in the file my_blast.xml, all you need to do is to open the file for reading:

# Save result into my_blast.xml and read that file 
# result_handle = open("my_blast.xml")

# ================================================================================
# @@ Now that we’ve got a handle, we are ready to parse the output. 

# ================================================================================
# The code to parse it is really quite small. If you expect a single BLAST result (i.e., you used a single query):

# result_handle: handle containing my_blast.xml
# blast_record = NCBIXML.read(result_handle)

# ================================================================================
# or, if you have lots of results (i.e., multiple query sequences):

# Parse result_handle
# blast_records = NCBIXML.parse(result_handle)

# ================================================================================
# Just like Bio.SeqIO and Bio.AlignIO (see Chapters 5 and 6), we have a pair of input functions, read and parse, where read is for when you have exactly one object, and parse is an iterator for when you can have lots of objects – but instead of getting SeqRecord or MultipleSeqAlignment objects, we get BLAST record objects.

# To be able to handle the situation where the BLAST file may be huge, containing thousands of results, NCBIXML.parse() returns an iterator. 

# ================================================================================
# In plain English, an iterator allows you to step through the BLAST output, retrieving BLAST records one by one for each BLAST search result:

# ================================================================================
# c blast_records: iterator which contains parsed data from result_handle
# blast_records = NCBIXML.parse(result_handle)

# blast_record = next(blast_records)
# ... do something with blast_record

# blast_record = next(blast_records)
# ... do something with blast_record

# blast_record = next(blast_records)
# ... do something with blast_record

# blast_record = next(blast_records)
# Traceback (most recent call last):
#   File "<stdin>", line 1, in <module>
# StopIteration
# No further records

# ================================================================================
# Or, you can use a for-loop:

# for blast_record in blast_records:
#     Do something with blast_record

# ================================================================================
# Note though that you can step through the BLAST records only once. Usually, from each BLAST record you would save the information that you are interested in. 

# ================================================================================
# If you want to save all returned BLAST records, you can convert the iterator into a list:

# c blast_records: "blast_records iterator" into list
# blast_records = list(blast_records)

# ================================================================================
# @@ Now you can access each BLAST record in the list with an index as usual. 

# ================================================================================
# @@ If your BLAST file is huge though, you may run into memory problems trying to save them all in a list.

# ================================================================================
# Usually, you’ll be running one BLAST search at a time. Then, all you need to do is to pick up the first (and only) BLAST record in blast_records:

# blast_records = NCBIXML.parse(result_handle)
# blast_record = next(blast_records)

# ================================================================================
# or more elegantly:
# blast_record = NCBIXML.read(result_handle)

# ================================================================================
# I guess by now you’re wondering what is in a BLAST record.
