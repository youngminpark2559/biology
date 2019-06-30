# http://www.incodom.kr/Biopython/Bio.AlignIO

# 계통 분석 등을 하기 위해 자주 사용하는 툴은 ClustalW가 있는데 이에 대한 결과를 파싱하고 활용하기 위해 직접 데이터를 핸들링하기에는 어려움이 따른다.

# 더구나 방대한 데이터를 메뉴얼로 분석한다는 것은 상상하기 힘들다.

# 그래서 biopython에서는 AlignIO 패키지를 제공하고 있다.

# 이는 다양한 포맷의 정렬 결과를 읽고 데이터를 추출하는 기능을 제공하고 있다.

# ================================================================================
# 먼저, 사용을 위해 from Bio import AlignIO와 같이 모듈을 임포트 하면 된다.

# 그 후 사용법은 help(모듈명)를 통해 확인해 볼 수 있는데, 도움말을 통해 어떤 기능이 있는지 확인 후에 이용하면 된다.

from Bio import AlignIO
help(AlignIO)

# ================================================================================
# Bio.SeqIO 패키지와 같이 .read()와 .parse()를 사용할 수 있다.

# 이는 단일 결과와 다중 결과를 처리한다는 차이점이 존재한다.

# 정렬 포맷에는 다양한 형태가 있는데 그 중 첫번째로 스톡홀름 포맷은 다음과 같이 구성되어 있다.

# 입력한 서열의 정보와 정렬된 결과 서열 및 위치 정보를 출력하는 두 부분으로 나누어진다.

# 다음 절에서 이 포맷의 파일을 읽어 들이는 예제를 확인할 수 있다.

# # STOCKHOLM 1.0
# #=GS COATB_BPIKE/30-81  AC P03620.1
# #=GS COATB_BPIKE/30-81  DR PDB; 1ifl ; 1-52;
# #=GS Q9T0Q8_BPIKE/1-52  AC Q9T0Q8.1
# #=GS COATB_BPI22/32-83  AC P15416.1
# #=GS COATB_BPM13/24-72  AC P69541.1
# #=GS COATB_BPM13/24-72  DR PDB; 2cpb ; 1-49;
# #=GS COATB_BPM13/24-72  DR PDB; 2cps ; 1-49;
# #=GS COATB_BPZJ2/1-49   AC P03618.1
# #=GS Q9T0Q9_BPFD/1-49   AC Q9T0Q9.1
# #=GS Q9T0Q9_BPFD/1-49   DR PDB; 1nh4 A; 1-49;
# #=GS COATB_BPIF1/22-73  AC P03619.2
# #=GS COATB_BPIF1/22-73  DR PDB; 1ifk ; 1-50;
# COATB_BPIKE/30-81             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
# #=GR COATB_BPIKE/30-81  SS    -HHHHHHHHHHHHHH--HHHHHHHH--HHHHHHHHHHHHHHHHHHHHH----
# Q9T0Q8_BPIKE/1-52             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
# COATB_BPI22/32-83             DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
# COATB_BPM13/24-72             AEGDDP...AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
# #=GR COATB_BPM13/24-72  SS    ---S-T...CHCHHHHCCCCTCCCTTCHHHHHHHHHHHHHHHHHHHHCTT--
# COATB_BPZJ2/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
# Q9T0Q9_BPFD/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
# #=GR Q9T0Q9_BPFD/1-49   SS    ------...-HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--
# COATB_BPIF1/22-73             FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA
# #=GR COATB_BPIF1/22-73  SS    XX-HHHH--HHHHHH--HHHHHHH--HHHHHHHHHHHHHHHHHHHHHHH---
# #=GC SS_cons                  XHHHHHHHHHHHHHHHCHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHHC--
# #=GC seq_cons                 AEssss...AptAhDSLpspAT-hIu.sWshVsslVsAsluIKLFKKFsSKA
# //

# ================================================================================
# @@ 모듈 구성을 불러오기 위해서는 다음과 같이 호출한다.

# @@ 호출은 기존과 마찬가지로 Bio 패키지에서 AlignIO 모듈을 불러 오면 된다.

# @@ 모듈을 임포트 후에 다음과 같은 방법으로 파일을 읽어 들일 수 있다.

# 여기서는 위의 파일 포맷과 같이 스톡홀름 포맷(Stockholm format)이므로 아래와 같이 지정한다.

# 즉, AlignIO.read("파일이름", "stockholm")와 같이 사용하면 된다.

from Bio import AlignIO
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

# ================================================================================
# @@ 스톡홀름 포맷으로 파일을 읽었으면 아래와 같이 출력해 볼 수 있다. 

# 파일을 읽어서 표준 alignment에 대한 형태로 출력을 해준다.
print(alignment)

# SingleLetterAlphabet() alignment with 7 rows and 52 columns
# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
# DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
# AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
# FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

# ================================================================================
# 마찬가지 방법으로 각각 alignment 된 전체 결과 사이즈를 구할 수 있다.

# 이는 alignment.get_alignment_length() 함수를 통해 구할 수 있다.

# 더불어 순환문을 통해 각 alignment의 서열 및 서열명을 추출해 낼 수 있다.

# 아래 예제에서는 그 방법을 설명하고 있다.

from Bio import AlignIO
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
# print("Alignment length %i" % alignment.get_alignment_length())
# Alignment length 52

for record in alignment:
    print("%s - %s" % (record.seq, record.id))

# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
# DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
# AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
# FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

# ================================================================================
# 또한 다음과 같이 FASTA 형태의 결과를 가지고 있을 수도 있다.

# >COATB_BPIKE/30-81
# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
# >Q9T0Q8_BPIKE/1-52
# AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
# >COATB_BPI22/32-83
# DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
# >COATB_BPM13/24-72
# AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
# >COATB_BPZJ2/1-49
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
# >Q9T0Q9_BPFD/1-49
# AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
# >COATB_BPIF1/22-73
# FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA

# ================================================================================
# 이 파일은 다음과 같이 FASTA 형태로 데이터를 읽으면 파일을 읽을 수가 있다. 

# AlignIO.read("PF05371_seed.faa", "fasta")를 이용하여 포맷형태만 바꾸어 주면 된다.

from Bio import AlignIO
alignment = AlignIO.read("PF05371_seed.faa", "fasta")
print(alignment)

# ================================================================================
# 다음의 예제는 phylip 포맷이다. 

# 역시 계통분석을 위해서 사용하는 phylip이나 clustalw에서 phylip 형태로 저장할 수 있다. 

# 아래 포맷도 마찬가지로 서열명, 서열 정보를 가지고 있다.

#     5    6
# Alpha     AACAAC
# Beta      AACCCC
# Gamma     ACCAAC
# Delta     CCACCA
# Epsilon   CCAAAC

# ================================================================================
# 이 경우에는 역시 AlignIO.parse("resampled.phy", "phylip")와 같이 파일 이름과 파일 포맷 형태만 바꿔주면 된다. 

# 나머지 사용법은 동일하다. 

# 역시 순환문을 통해 각 정렬된 서열을 출력할 수 있다.

from Bio import AlignIO
alignments = AlignIO.parse("resampled.phy", "phylip")
for alignment in alignments:
    print(alignment)
    print("")

# ================================================================================
# 위에서 다루었던 부분은 정렬 파일을 읽어 활용을 했다면 이제부터는 정렬 결과를 파일에 저장하는 예제이다.

# 물론 읽은 정렬 결과를 바로 저장할 수도 있지만 아래의 예제에서는 임의의 정렬 결과를 만드는 예제이다.

# 먼저 서열을 다루기 위해서 from Bio.Alphabet import generic_dna 모듈과 from Bio.Seq import Seq를 임포트 한다.

# 파일에 쓰기 위해 바이오파이썬에서 제공하는 서열레코드(SeqRecord)기능을 이용하면 된다.

# 그리고 마지막으로 MultipleSeqAlignment를 임포트 한다.

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# @@ 다음에는 각 정렬 파일을 만드는 과정이다. 

# @@ 정렬 파일은 MultipleSeqAlignment를 이용하여 리스트를 만들어 주면 된다. 

# @@ 각 하나의 서열을 만들고 그 서열을 레코드 객체로 만들어 처리하면 된다. 

# @@ 사용법은 다음과 같다. 

# SeqRecord(Seq("서열", 서열타입), id="서열명") 

# 만들어진 서열레코드를 정렬파일화 시키는 작업을 하고 그 정렬 결과를 파일에 쓰면 된다. 

# 주의할 점은 각 정렬 파일은 이중 리스트 구조를 가지고 있다는 점을 주의해야 한다.

align1 = MultipleSeqAlignment(
    [SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
     SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
     SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),])

align2 = MultipleSeqAlignment(
    [SeqRecord(Seq("GTCAGC-AG", generic_dna), id="Delta"),
     SeqRecord(Seq("GACAGCTAG", generic_dna), id="Epsilon"),
     SeqRecord(Seq("GTCAGCTAG", generic_dna), id="Zeta"),])

align3 = MultipleSeqAlignment(
    [SeqRecord(Seq("ACTAGTACAGCTG", generic_dna), id="Eta"),
     SeqRecord(Seq("ACTAGTACAGCT-", generic_dna), id="Theta"),
     SeqRecord(Seq("-CTACTACAGGTG", generic_dna), id="Iota"),])

my_alignments = [align1, align2, align3]

# @@ 이렇게 만들어진 정렬 결과는 SeqIO와 같이 동일한 방법으로 사용하면 된다. 

# 사용법은 AlignIO.write(정렬레코드, "파일이름", "파일포맷")이다. 

# 여기서는 phylip 형태로 저장하기로 한다.

from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")

# ================================================================================
# @@ 우리가 다루는 정렬 파일은 다양한 형태로 존재하고 필요에 따라 데이터의 포맷을 전환할 필요가 있다.

# @@ 이를 위해 간단히 전환하는 방법을 소개한다.

# @@ 사용법은 다음과 같다.

# AlignIO.convert("입력파일", "입력포맷", "출력파일", "출력포맷")와 사용을 하면 된다.

# 아래의 예제에서는 스톡홀롬 포맷에서 clustal 포맷으로 전환하는 예제이다.

from Bio import AlignIO
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
print("Converted %i alignments" % count)
