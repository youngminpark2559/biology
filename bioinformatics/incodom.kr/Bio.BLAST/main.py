# BLAST

# 바이오파이썬을 이용하여 NCBI에 BLAST를 수행할 수 있다. 

# 웹에서 이용하는 환경과 동일하게 이용할 수 있다. 

# 물론 직접 웹 브라우저를 이용하여 서열 분석을 할 수 도 있지만 갯수가 많거나 하면 많이 번거로울 수 밖에 없다. 

# 이는 바이오파이썬을 이용하면 간단하게 해결 할 수 있다.

# ================================================================================
from Bio.Blast import NCBIWWW
help(NCBIWWW.qblast)

# @@ SeqIO, AlignIO와 마찬가지로 단순히 패키지를 import 하면 된다.

# @@ from Bio.Blast import NCBIWWW 구문을 이용하여 사용한다.

# ================================================================================
# 아래의 예제에서는 nt 데이터베이스에 blastn을 하는데 서열은 "8332116" 이라는 GI를 가진 서열이다.

# 즉 NCBI에서 해당 서열을 찾아서 nt 데이터베이스에 검색을 하라는 것이다.

from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

# ================================================================================
# 물론 아래의 경우처럼 GI대신에 해당하는 서열을 입력해줘도 된다.

# 우리가 주로 사용하게 되는 서열을 입력하게 되면 자동으로 검색을 해서 결과를 알려준다.

from Bio.Blast import NCBIWWW
fasta_string = open("m_cold.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

# ================================================================================
# 아래의 예제에서는 해당 SeqIO패키지를 이용하여 FASTA 파일로 된 서열을 읽어서 해당 쿼리문을 만들 수 있다.

from Bio.Blast import NCBIWWW
from Bio import SeqIO

record = SeqIO.read("m_cold.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

# ================================================================================
# @@ 아래의 예제에서는 해당 내용을 파일로 저장하기 위한 방법이다.

# @@ 먼저 해당 내용을 저장할 파일의 핸들을 만들고 result_handle의 내용을 읽어서 저장하면 끝난다.

# 즉, 저장 파일 핸들 열기, BLAST 핸들 읽기, 그리고, 파일을 쓰는 과정을 거치면 된다.

save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()

# ================================================================================
# 조금 다르게 SeqIO의 기능을 이용하여 .format(포맷)기능을 이용하여 서열을 입력해 줄 수 있으며, 이는 원하는 형태로 사용하기 쉬운 방법이나 익숙한 방법으로 처리하면 된다.

# NCBIWWW의 기능을 이용하여 서열을 검색하고 비교 할 수 있지만 많은 경우 BLAST+가 설치된 로컬 서버에서도 사용한다.

# 즉 blastn의 명령어를 사용하여 직접 검색을 하고 결과를 가지고 있는 경우에는 커맨드 라인을 이용할 수 있다.

# 이 경우에는 from Bio.Blast.Applications import NcbiblastxCommandline 를 이용하여 해당 팩키지를 불러온다.

# 각 BLAST프로그램 용도에 맞게 구분되어 있으니 해당 용도에 맞게 임포트 해야 한다.

# 사용법은 커맨드 라인으로 사용할 때의 옵션과 별반 다르지 않다.

# blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db='nt') 등과 같이 사용하면 된다.

# 먼저 인스턴스를 하나 만들고 커맨드라인 NcbiblastxCommandline()에 파라미터를 입력해주면 된다.

from Bio.Blast.Applications import NcbiblastxCommandline
help(NcbiblastxCommandline)

blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")
# blastx_cline

NcbiblastxCommandline(cmd='blastx', out='opuntia.xml', outfmt=5, query='opuntia.fasta', db='nr', evalue=0.001)
print(blastx_cline)
# blastx -out opuntia.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001

stdout, stderr = blastx_cline()

# ================================================================================
그 이후 print(blastx_cline) 명령어를 통해 실제 로컬로 실행되어질 명령어를 확인할 수 있다.

이 명령어를 직접 실행해도 되지만 바이오파이썬에서는 stdout, stderr = blastx_cline() 와 같이 간단히 실행할 수 있다.

이는 표준 출력과 에러 출력으로 구분하여 결과를 저장해주며 필요에 따라 저장하면 된다.

BLAST를 사용하면서 그 결과는 넷 핸들과 파일 핸들로 나눌 수 있다.

NCBIWWW와 NCBIXML가 그것인데 전자가 넷 핸들이며 후자가 파일 핸들이다.

전자는 직접 브라우저를 호출하는 것이며 후자는 결과 파일 특히 XML 형식으로 저장된 데이터를 읽어서 결과를 파싱할 수 있다.

from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)

# 물론 여기에서도 read()와 parse() 의 기능은 SeqIO 와 마찬가지로 단일 혹은 멀티 결과 유무에 따라 호출하는 함수도 달라져야 한다.

# 아래의 예제에서는 HSP 영역을 추출해 내는 스크립트이다.

# 위의 결과에서 blast_records ==NCBIXML.parse(result_handle)로 호출하였을 때에는 해당 내용을 이용하기 위해서 loop 구문을 한번 더 돌아야 한다.

# 즉 for blast_record in blast_records : 와 같이 문이 상위에 들어가야 한다.

# 각 해당 내용을 확인하기 위해서는 dir(blast_record), dir(alignment), dir(hsp) 와 같이 해당 객체에 어떤 값 및 함수가 있는지 확인 한 후 해당 내용을 출력할 수 있다.

E_VALUE_THRESH = 0.04
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:75] + '...')
            print(hsp.match[0:75] + '...')
            print(hsp.sbjct[0:75] + '...')
