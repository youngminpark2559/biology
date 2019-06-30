# conda activate py36gputorch100 && \
# cd /mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence && \
# rm e.l && python main.py \
# 2>&1 | tee -a e.l && code e.l

# ================================================================================
# Original source:
# https://partrita.github.io/posts/python_randomseq/

# ================================================================================
from random import Random
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
Entrez.email="youngminpark2559@gmail.com"

# ================================================================================
def save_blast_result_into_xml(blast_result_xml_file_path):
  with open(blast_result_xml_file_path, "w") as out_handle:
    out_handle.write(result_handle.read())
  result_handle.close()

def load_blast_result_from_xml(blast_result_xml_file_path):
  result_handle = open(blast_result_xml_file_path)
  # print("result_handle",result_handle)
  return result_handle

def see_help_text_of_NCBIWWW_qblast():

  help(NCBIWWW.qblast)

  # Help on function qblast in module Bio.Blast.NCBIWWW:

  # qblast(program, database, sequence, url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', auto_format=None, composition_based_statistics=None, db_genetic_code=None, endpoints=None, entrez_query='(none)', expect=10.0, filter=None, gapcosts=None, genetic_code=None, hitlist_size=50, i_thresh=None, layout=None, lcase_mask=None, matrix_name=None, nucl_penalty=None, nucl_reward=None, other_advanced=None, perc_ident=None, phi_pattern=None, query_file=None, query_believe_defline=None, query_from=None, query_to=None, searchsp_eff=None, service=None, threshold=None, ungapped_alignment=None, word_size=None, alignments=500, alignment_view=None, descriptions=500, entrez_links_new_window=None, expect_low=None, expect_high=None, format_entrez_query=None, format_object=None, format_type='XML', ncbi_gi=None, results_file=None, show_overview=None, megablast=None, template_type=None, template_length=None)
  #     BLAST search using NCBI's QBLAST server or a cloud service provider.
      
  #     Supports all parameters of the qblast API for Put and Get.
      
  #     Please note that BLAST on the cloud supports the NCBI-BLAST Common
  #     URL API (http://ncbi.github.io/blast-cloud/dev/api.html). To
  #     use this feature, please set url_base to
  #     'http://host.my.cloud.service.provider.com/cgi-bin/blast.cgi' and
  #     format_object='Alignment'. For more details, please see
  #     https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=CloudBlast
      
  #     Some useful parameters:
      
  #      - program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
  #      - database       Which database to search against (e.g. "nr").
  #      - sequence       The sequence to search.
  #      - ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
  #      - descriptions   Number of descriptions to show.  Def 500.
  #      - alignments     Number of alignments to show.  Def 500.
  #      - expect         An expect value cutoff.  Def 10.0.
  #      - matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
  #      - filter         "none" turns off filtering.  Default no filtering
  #      - format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
  #      - entrez_query   Entrez query to limit Blast search
  #      - hitlist_size   Number of hits to return. Default 50
  #      - megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
  #      - service        plain, psi, phi, rpsblast, megablast (lower case)
      
  #     This function does no checking of the validity of the parameters
  #     and passes the values to the server as is.  More help is available at:
  #     https://ncbi.github.io/blast-cloud/dev/api.html

  # nt	Nucleotide collection	DNA
  # nr	Non-redundant	Protein
  # refseq_rna	NCBI Transcript Reference Sequences	DNA
  # refseq_protein	NCBI Protein Reference Sequences	Protein
  # swissprot	Non-redundant UniProtKB/SwissProt sequences	Protein
  # pdbaa	PDB protein database	Protein
  # pdbnt	PDB nucleotide database	DNA

def BLAST_by_using_NCBIWWW_qblast(blast_program,database_for_blast,your_data):

  # record = SeqIO.read("m_cold.fasta", format="fasta")
  result_handle = NCBIWWW.qblast(blast_program,database_for_blast,your_data)

  # print("result_handle",result_handle)
  # <_io.StringIO object at 0x7fd5b7ee2d38>

  return result_handle

def AminoacidGenerator(number, length, output_text_file_path):
    print('***generate random aminoacid sequence ***')

    # ================================================================================
    codeFile = open(output_text_file_path, 'w')
    if number <= 0:
        return 'invalid number of protein'
    else:
        aminoacid = 'ACDEFGHIKLMNPQRSTVWY'      
        random = Random()


        all_random_seqs=[]
        for j in range(1, number + 1):
            random_aa_str = ''
            for i in range(1, length+1):
                index = random.randint(1,len(aminoacid))
                random_aa_str = random_aa_str + aminoacid[index -1]
                # print("random_aa_str",random_aa_str)
                # R

            # print("random_aa_str",random_aa_str)
            # CNHTRFADFAGGIPLKWMNSMYRAPCDENCFSYQGYGSEWDHIFLLQDYMAGDTKILWGWREMSEGVGKQYVNKTCSHVEMLYYPQCNEPHLADIPLNWT

            all_random_seqs.append(random_aa_str)

            codeFile.write(random_aa_str+'\n')

    return all_random_seqs

# ================================================================================
# c num_aa_seq: number of amino acide sequence
num_aa_seq=100

# c len_aa_seq: length of amino acid sequence
len_aa_seq=100

output_text_file_path="/mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/generated_random_amino_acid_sequence_"+str(len_aa_seq)+"_"+str(num_aa_seq)+".txt"
# print("output_text_file_path",output_text_file_path)
# /mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/generated_random_amino_acid_sequence_100_100.txt

blast_result_xml_file_path="/mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/my_blast.xml"

# ================================================================================
all_random_seqs=AminoacidGenerator(num_aa_seq,len_aa_seq,output_text_file_path)
# print(all_random_seqs)
# ['FVHYMKAQDPFAHPTDVTHPFNYKMQPCYSESNWKWQKAWMWNFSFGPRDRNWYLYAPKFVFRFSWNKNAWMMMWIVAGEESNNEPQPIHKIASKQPKHM', 
#  'PSLHEWPMYANAHESQSEMWKEKIFYNCNGGEAQMFNYIDDACGAVWCKNMQFFLTWTVMQVHGMNIFFRVWQQWMGFCCKYIYGTTLDNYQKNCYQDFR', 
#  'YEPLDDVTIPSCNMNNEYSQVPNKNMVKQQWYLFFACSKLNECRIQRNGLRLKQDVEMQLSEIAHKEDKRFECMHVLDWFGSMGWEGAHGTQFDSFPIAW', 

# ================================================================================
# @ See the help text about NCBIWWW.qblast() if you need it
# see_help_text_of_NCBIWWW_qblast()

# ================================================================================
# @ How to do BLAST using Web BLAST
# - <https://blast.ncbi.nlm.nih.gov/Blast.cgi>
# - Click Protein BLAST
# - You're redirected to
# - <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome>
# - Paste following sequence
# MSKRKAPQETLNGGITDMLTELANFEKNVSQAIHKYNAYRKAASVIAKYPHKIKSGAEAKKLPGVGTKIAEKIDEFLATGKLRKLEKIRQDDTSSSINFLTRVSGIGPSAARKFVDEGIKTLEDLRKNEDKLNHHQRIGLKYFGDFEKRIPREEMLQMQDIVLNEVKKVDSEYIATVCGSFRRGAESSGDMDVLLTHPSFTSESTKQPKLLHQVVEQLQKVHFITDTLSKGETKFMGVCQLPSKNDEKEYPHRRIDIRLIPKDQYYCGVLYFTGSDIFNKNMRAHALEKGFTINEYTIRPLGVTGVAGEPLPVDSEKDIFDYIQWKYREPKDRSE
# - Click BLAST

# ================================================================================
# one_aa="MSKRKAPQETLNGGITDMLTELANFEKNVSQAIHKYNAYRKAASVIAKYPHKIKSGAEAKKLPGVGTKIAEKIDEFLATGKLRKLEKIRQDDTSSSINFLTRVSGIGPSAARKFVDEGIKTLEDLRKNEDKLNHHQRIGLKYFGDFEKRIPREEMLQMQDIVLNEVKKVDSEYIATVCGSFRRGAESSGDMDVLLTHPSFTSESTKQPKLLHQVVEQLQKVHFITDTLSKGETKFMGVCQLPSKNDEKEYPHRRIDIRLIPKDQYYCGVLYFTGSDIFNKNMRAHALEKGFTINEYTIRPLGVTGVAGEPLPVDSEKDIFDYIQWKYREPKDRSE"
# my_prot=Seq(one_aa,IUPAC.protein)
# blast_result_try_00001=BLAST_by_using_NCBIWWW_qblast(blast_program="blastp",database_for_blast="nr",your_data=my_prot)
# print("blast_result_try_00001",blast_result_try_00001.read())
# /mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/blast_result_try_00001.xml

# ================================================================================
# one_aa="MDSTGEFCWICHQPEGPLKRFCGCKGSCAVSHQDCLRGWLETSRRQTCALCGTPYSMKWKTKPLREWTWGEEEVLAAMEACLPLVLIPLAVLMIVMGTWLLVNHNGFLSPRMQVVLVVIVLLAMIVFSASASYVMVEGPGCLDTCTAKNSTVTVNSIDEAIATQQPTKTDLGLARETLSTRFRRGKCRSCCRLGCVRLCCV"
# my_prot=Seq(one_aa,IUPAC.protein)
# blast_result_try_00002=BLAST_by_using_NCBIWWW_qblast(blast_program="blastp",database_for_blast="nr",your_data=my_prot)
# print("blast_result_try_00002",blast_result_try_00002.read())
# /mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/blast_result_try_00002.xml

# ================================================================================
# one_aa=all_random_seqs[0]
# my_prot=Seq(one_aa,IUPAC.protein)
# blast_result_try_00003=BLAST_by_using_NCBIWWW_qblast(blast_program="blastp",database_for_blast="nr",your_data=my_prot)
# print("blast_result_try_00003",blast_result_try_00003.read())
# /mnt/1T-5e7/mycodehtml/Biology/partrita.github.io/Create_random_protein_sequence/Results_log/blast_result_try_00003.xml

# ================================================================================
# save_blast_result_into_xml(blast_result_xml_file_path):

# ================================================================================
# load_blast_result_from_xml(blast_result_xml_file_path):
