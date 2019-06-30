# http://www.incodom.kr/Biopython/Bio.Entrez

# 바이오파이썬을 이용하여 NCBI데이터 핸들링하기 #
# 미국 국립생물공학정보센터(National Center for Biotechnology Information(NCBI))는 미국 보건성 산하 국립의학도서관의 운영 분야 중 하나이다. NCBI는 1988년 미국 메릴랜드 주에 설립되었으며 생명과학 및 의학 논문 인덱스의 데이터베이스인 펍메드(PubMed), 유전체 서열 데이터베이스인 진뱅크(GenBank)를 비롯하여 각종 생명공학 정보들을 담고 있으며, 이 모든 정보들은 Entrez 검색엔진을 이용하여 온라인으로 열람할 수 있다. 이러한 검색엔진인 Entrez를 이용할 수 있도록 각종 Bio**** 프로그램들이 이 모듈을 탑재하고 있다. 바이오파이썬에도 마찬가지로 Entrez를 이용하여 데이터 검색 및 다운로드할 수 있다. 그중에서도 주요 데이터의 종류 검색(einfo), 데이터 검색(esearch), 검색한 결과를 확인 및 저장하는 방법(efetch)에 대해 알아보고자 한다.

from Bio import Entrez
# 먼저 from Bio import Entrez 팩키지를 이용한다. 이렇게 불러 들이면 준비가 된 것이다. 주의할 사항은 Entrez.email='test@example.org' 와 같이 메일 주소를 적어주지 않으면 에러 메시지를 뿌려준다.

# 아래와 같이 반드시 메일 주소를 알려주어야 한다. einfo를 이용하여 Entrez에서 제공하는 데이터베이스 정보를 확인할 수 있다. 사용법은 Entrez.einfo()를 호출하면 되며 이를 핸들 정보에 저장하고 핸들에 저장된 정보를 handle.read()를 이용하면 된다.

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"

# ================================================================================
# EINFO

# @@ 핸들에 저장된 정보를 출력하면 다음과 같이 출력된다.

# @@ 기본적으로 XML(데이터 교환 포맷) 형태로 출력하게 되는데 이 안에 데이터베이스의 종류를 확인할 수 있다.

# @@ 흔히 알고 있는 pubmed부터 각종 서열 정보까지 NCBI에서 제공하는 모든 리스트를 확인할 수 있다.

handle = Entrez.einfo()
result = handle.read()
print(result)
# <?xml version="1.0"?>
# <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN"
#  "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
# <eInfoResult>
# <DbList>
#         <DbName>pubmed</DbName>
#         <DbName>protein</DbName>
#         <DbName>nucleotide</DbName>
#         <DbName>nuccore</DbName>
#         <DbName>nucgss</DbName>
#         <DbName>nucest</DbName>
#         <DbName>structure</DbName>
#         <DbName>genome</DbName>
#         <DbName>books</DbName>
#         <DbName>cancerchromosomes</DbName>
#         <DbName>cdd</DbName>
#         <DbName>gap</DbName>
#         <DbName>domains</DbName>
#         <DbName>gene</DbName>
#         <DbName>genomeprj</DbName>
#         <DbName>gensat</DbName>
#         <DbName>geo</DbName>
#         <DbName>gds</DbName>
#         <DbName>homologene</DbName>
#         <DbName>journals</DbName>
#         <DbName>mesh</DbName>
#         <DbName>ncbisearch</DbName>
#         <DbName>nlmcatalog</DbName>
#         <DbName>omia</DbName>
#         <DbName>omim</DbName>
#         <DbName>pmc</DbName>
#         <DbName>popset</DbName>
#         <DbName>probe</DbName>
#         <DbName>proteinclusters</DbName>
#         <DbName>pcassay</DbName>
#         <DbName>pccompound</DbName>
#         <DbName>pcsubstance</DbName>
#         <DbName>snp</DbName>
#         <DbName>taxonomy</DbName>
#         <DbName>toolkit</DbName>
#         <DbName>unigene</DbName>
#         <DbName>unists</DbName>
# </DbList>
# </eInfoResult>

# ================================================================================
# @@ 이를 좀 더 편하게 확인하기 위해 바이오파이썬에는 DB목록을 파싱하여 제공한다.

# @@ 아래 예제처럼 읽어들인 핸들 정보를 레코드에 담아 순환문을 통해 정보를 활용할 수 있게 제공한다.

# @@ 파이썬의 자료 구조인 사전(dictionary)에 저장되어 있기 떼문에 먼저 record.keys()를 이용하여 키 목록을 확인하다.

# @@ 그 후 record['DbList'] 를 이용하여 값을 불러 온다.

# @@ 파이썬에서 사전 구조를 읽는 방법은 dict['key'] 의 형태로 사용하면 된다.

from Bio import Entrez
handle = Entrez.einfo()
record = Entrez.read(handle)
print(record.keys())
# [u'DbList']
print(record["DbList"])
# ['pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
#  'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
#  'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
#  'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
#  'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
#  'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists']

# ================================================================================
# 레코드에 저장된 정보 중 'DbInfo'에는 각 데이터에 대한 통계 정보를 담고 있다. 

# 그 중 논문 초록 정보를 조회하려면 'pubmed'를 이용하면 된다. 

# 즉 einfo(db="db_name")과 같이 사용하면 된다. 

# db에 위의 예제 record["DbList"]에서 확인한 데이터의 종류를 입력하면 된다.

handle = Entrez.einfo(db="pubmed")
record = Entrez.read(handle)
record["DbInfo"]["Description"]
# 'PubMed bibliographic record'
record["DbInfo"]["Count"]
# '17989604'
record["DbInfo"]["LastUpdate"]
# '2008/05/24 06:45

# ================================================================================
# 그러면 이제 pubmed 데이터베이스에 대한 정보를 확인할 수 있다.

# 설명(Description), 데이터 수(Count), 최근 업데이트 날짜(LastUpdate)와 같은 정보를 확인하면 된다.

# 위의 예제에서 그 사용법을 보여 주고 있다.

# ================================================================================
# ESEARCH

# 다음으로는 방대한 자료를 검색할 수 있는 방법을 알아보고자 한다.

# 검색은 아래 예제에서처럼 esearch()를 이용한다.

# 사용법은 Entrez.esearch(db="데이터베이스이름", term="검색어")와 같이 사용한다.

# db에는 'DbList'에서 찾은 데이터베이스 이름을 입력하고 term에는 우리가 찾아야 할 키워드를 입력하면 된다.

# 아래 예제에서는 pubmed에서 biopython이라는 키워드가 들어간 논문 정보를 찾겠다는 내용이다.

# 마찬가지로 정보를 핸들(handle) 변수에 저장하고 정보를 읽어(read(handle)) 레코드에 저장한다.

# 레코드 변수에는 역시 사전 구조로 데이터를 저장하고 있다.

# 검색 결과의 고유 ID 및 찾은 갯수 등을 저장한다.

# 'IdList'를 통해 각 정보를 확인한다.

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
handle = Entrez.esearch(db="pubmed", term="biopython")
record = Entrez.read(handle)
record["IdList"]
# ['19304878', '18606172', '16403221', '16377612', '14871861', '14630660', '12230038']

# ================================================================================
# 다음 예제에서는 서열을 검색하기 위하여 db='nucleotide', 검색어에 term='Cypripedioideae[Orgn] AND matK[Gene]'를 넣었다.

# 검색어에는 정규표현식과 같은 복잡한 문법을 사용할 수 있다.

# 여기에서는 [Orgn]은 특정 계통(taxonomy)을 선택, [Gene]은 matK라는 유전자 이름을 가진 것만 찾겠다는 것이다.

# 역시 검색하고 정보를 읽어들이면 검색 결과를 확인할 수 있다.

handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]")
record = Entrez.read(handle)
record["Count"]
# '25'
record["IdList"]
# ['126789333', '37222967', '37222966', '37222965',    , '61585492']

# ================================================================================
# 다음 예제에서는 저널 정보를 검색하는 것이다.

handle = Entrez.esearch(db="journals", term="computational")
record = Entrez.read(handle)
record["Count"]
# '16'
record["IdList"]
# ['30367', '33843', '33823', '32989', '33190', '33009', '31986',
#  '34502', '8799', '22857', '32675', '20258', '33859', '32534',
#  '32357', '32249']

# ================================================================================
# EFETCH

# 다음으로는 검색한 결과를 세부 정보를 확인하는 방법이다.

# 사용법은 efetch()를 사용한다.

# 기본적으로 다음과 같이 efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text") 사용한다.

# esearch()와 마찬가지로 db(데이터 종류), 인덱스 번호(id),데이터 포맷(rettype), 데이터 형태(retmode)를 이용하면 된다.

# 여기에서는 특정 인덱스 번호 186972394(esearch에서 찾은 번호 혹은 이미 알고 있는 번호)를 서열데이터인 nucleotide, genbank포맷으로 읽겠다고 알려주는 것이다.

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"
handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
print(handle.read())
# LOCUS       EU490707                1302 bp    DNA     linear   PLN 05-MAY-2008
# DEFINITION  Selenipedium aequinoctiale maturase K (matK) gene, partial cds;
#             chloroplast.
# ACCESSION   EU490707
# VERSION     EU490707.1  GI:186972394
# KEYWORDS    .
# SOURCE      chloroplast Selenipedium aequinoctiale
#   ORGANISM  Selenipedium aequinoctiale
#             Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
#             Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae;
#             Cypripedioideae; Selenipedium.
# REFERENCE   1  (bases 1 to 1302)
#   AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
#             Endara,C.L., Williams,N.H. and Moore,M.J.
#   TITLE     Phylogenetic utility of ycf1 in orchids
#   JOURNAL   Unpublished
# REFERENCE   2  (bases 1 to 1302)
#   AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
#             Endara,C.L., Williams,N.H. and Moore,M.J.
#   TITLE     Direct Submission
#   JOURNAL   Submitted (14-FEB-2008) Department of Botany, University of
#             Florida, 220 Bartram Hall, Gainesville, FL 32611-8526, USA
# FEATURES             Location/Qualifiers
#      source          1..1302
#                      /organism="Selenipedium aequinoctiale"
#                      /organelle="plastid:chloroplast"
#                      /mol_type="genomic DNA"
#                      /specimen_voucher="FLAS:Blanco 2475"
#                      /db_xref="taxon:256374"
#      gene            <1..>1302
#                      /gene="matK"
#      CDS             <1..>1302
#                      /gene="matK"
#                      /codon_start=1
#                      /transl_table=11
#                      /product="maturase K"
#                      /protein_id="ACC99456.1"
#                      /db_xref="GI:186972395"
#                      /translation="IFYEPVEIFGYDNKSSLVLVKRLITRMYQQNFLISSVNDSNQKG
#                      FWGHKHFFSSHFSSQMVSEGFGVILEIPFSSQLVSSLEEKKIPKYQNLRSIHSIFPFL
#                      EDKFLHLNYVSDLLIPHPIHLEILVQILQCRIKDVPSLHLLRLLFHEYHNLNSLITSK
#                      KFIYAFSKRKKRFLWLLYNSYVYECEYLFQFLRKQSSYLRSTSSGVFLERTHLYVKIE
#                      HLLVVCCNSFQRILCFLKDPFMHYVRYQGKAILASKGTLILMKKWKFHLVNFWQSYFH
#                      FWSQPYRIHIKQLSNYSFSFLGYFSSVLENHLVVRNQMLENSFIINLLTKKFDTIAPV
#                      ISLIGSLSKAQFCTVLGHPISKPIWTDFSDSDILDRFCRICRNLCRYHSGSSKKQVLY
#                      RIKYILRLSCARTLARKHKSTVRTFMRRLGSGLLEEFFMEEE"
# ORIGIN
#         1 attttttacg aacctgtgga aatttttggt tatgacaata aatctagttt agtacttgtg
#        61 aaacgtttaa ttactcgaat gtatcaacag aattttttga tttcttcggt taatgattct
#       121 aaccaaaaag gattttgggg gcacaagcat tttttttctt ctcatttttc ttctcaaatg
#       181 gtatcagaag gttttggagt cattctggaa attccattct cgtcgcaatt agtatcttct
#       241 cttgaagaaa aaaaaatacc aaaatatcag aatttacgat ctattcattc aatatttccc
#       301 tttttagaag acaaattttt acatttgaat tatgtgtcag atctactaat accccatccc
#       361 atccatctgg aaatcttggt tcaaatcctt caatgccgga tcaaggatgt tccttctttg
#       421 catttattgc gattgctttt ccacgaatat cataatttga atagtctcat tacttcaaag
#       481 aaattcattt acgccttttc aaaaagaaag aaaagattcc tttggttact atataattct
#       541 tatgtatatg aatgcgaata tctattccag tttcttcgta aacagtcttc ttatttacga
#       601 tcaacatctt ctggagtctt tcttgagcga acacatttat atgtaaaaat agaacatctt
#       661 ctagtagtgt gttgtaattc ttttcagagg atcctatgct ttctcaagga tcctttcatg
#       721 cattatgttc gatatcaagg aaaagcaatt ctggcttcaa agggaactct tattctgatg
#       781 aagaaatgga aatttcatct tgtgaatttt tggcaatctt attttcactt ttggtctcaa
#       841 ccgtatagga ttcatataaa gcaattatcc aactattcct tctcttttct ggggtatttt
#       901 tcaagtgtac tagaaaatca tttggtagta agaaatcaaa tgctagagaa ttcatttata
#       961 ataaatcttc tgactaagaa attcgatacc atagccccag ttatttctct tattggatca
#      1021 ttgtcgaaag ctcaattttg tactgtattg ggtcatccta ttagtaaacc gatctggacc
#      1081 gatttctcgg attctgatat tcttgatcga ttttgccgga tatgtagaaa tctttgtcgt
#      1141 tatcacagcg gatcctcaaa aaaacaggtt ttgtatcgta taaaatatat acttcgactt
#      1201 tcgtgtgcta gaactttggc acggaaacat aaaagtacag tacgcacttt tatgcgaaga
#      1261 ttaggttcgg gattattaga agaattcttt atggaagaag aa
# //

# ================================================================================
# 바이오파이썬에서는 서열 데이터를 읽고 쓸 수 있다.

# SeqIO 팩키지를 이용하면 되는데 먼저 efetch()를 통해 검색하고 여기에서 읽은 정보를 SeqIO의 입력으로 주면 된다.

# 검색 결과 포맷은 gb(genbank)로 가져왔기 때문에 SeqIO()에서도 마찬가지로 파일 포맷 형태를 genbank라고 알려주면 된다.

# 이렇게 되면 위에 결과를 쉽게 핸들링(파싱) 할 수 있는 준비가 되어 보다 편하게 이용할 수 있다.

# 아래 예제에서는 검색하고 검색결과 서열을 핸들링하는 예제이다.

from Bio import Entrez, SeqIO
handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
print(record)
# ID: EU490707.1
# Name: EU490707
# Description: Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast.
# Number of features: 3
#    
# Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA   GAA', IUPACAmbiguousDNA())

# ================================================================================
# 다음의 예제에서는 검색 결과를 파일로 저장하는 방법을 보여주고 있다.

# 위의 에제에서는 임시로 프로그램 메모리에 저장되어 프로그램이 종료되면 결과를 확인할 수 없기 때문에 하드 디스크에 보관하기 위해서 검색 정보를 파일로 저장하는 방법을 보여주고 있다.

# 아래 예제에서는 간략하게 os 팩캐지를 이용해서 if not os.path.isfile(filename) 처럼 내가 저장하고자 하는 경로에 파일이 없으면 저장하라는 것이다.

# 그러면 파일 핸들을 열어 검색 결과를 저장하면 된다.

# 순환문을 이용하면 검색한 많은 정보를 자동으로 손쉽게 저장할 수 있다.

import os
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "A.N.Other@example.com"
filename = "gi_186972394.gbk"
if not os.path.isfile(filename):        # Downloading   
    net_handle = Entrez.efetch(db="nucleotide",id="186972394",rettype="gb", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()

print("Saved")
# Saved
print("Parsing   ")
# Parsing   
record = SeqIO.read(filename, "genbank")
print(record)
# ID: EU490707.1
# Name: EU490707
# Description: Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast.
# Number of features: 3
# /source=chloroplast Selenipedium aequinoctiale
# /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta',     'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Selenipedium']
# /keywords=['']
# /references=[Reference(title='Phylogenetic utility of ycf1 in orchids: a plastid gene more variable than matK',    ), Reference(title='Direct Submission',    )]
# /accessions=['EU490707']
# /data_file_division=PLN
# /date=26-JUL-2016
# /organism=Selenipedium aequinoctiale
# /sequence_version=1
# /topology=linear
# Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA   GAA', IUPACAmbiguousDNA())

