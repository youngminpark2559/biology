# Biopython을 이용한 Seq 객체 다루기

# @@ 우리가 컴퓨터로 데이터 특히 생물정보 서열을 다루기 위해서는 많은 노력이 필요하다.

# @@ 물론 간략하게 텍스트 에디터를 이용할 수 도 있고 GUI 프로그램을 이용할 수도 있다.

# @@ 그러나 이러한 어플리케이션에는 분명 한계가 존재하기 때문에 python 프로그래밍을 이용하여 해결할 수도 있다.

# @@ 이를 이용하기 위해서는 많은 문법 및 함수 등을 익히고 있어야 하므로 초보자가 사용하기에는 쉽지 않다.

# @@ 이때 Biopython의 서열을 다루기 위한 패키지 모듈을 이용하면 간단하게 해결할 수 있다.

# @@ 여기서는 python idle을 이용하여 진행한다.

# @@ 사용법은 다음과 같다.

# @@ from Bio.Seq import Seq 구문만 프로그램의 시작 부분에 선언하면 된다.

from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print("my_seq",my_seq)
# Seq('AGTACACTGGT', Alphabet())

print("my_seq.alphabet",my_seq.alphabet)
# Alphabet()

# ================================================================================
# @@ 처음에 my_seq라는 변수를 선언하고 Seq 객체를 이용하여 인스턴스를 생성하면 Seq 패키지 안에 선언된 다양한 함수 및 기능을 이용할 수 있다.

# @@ 사용할 수 있는 함수는 dir(my_seq) 를 이용하여 확인할 수 있다.

# @@ idle은 대화형으로 진행되기 때문에 바로바로 결과를 확인할 수 있다는 장점이 있다.

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
print("my_seq",my_seq)
# Seq('AGTACACTGGT', IUPACUnambiguousDNA())

print("my_seq.alphabet",my_seq.alphabet)
# IUPACUnambiguousDNA()

# ================================================================================
# 위의 예제에서는 from Bio.Alphabet import IUPAC을 이용하여 서열의 타입을 정의하여 필요한 사항을 사전 정의할 수 있다. 

# 사용법에는 Seq(서열, 타입코드)를 선언하면 된다. 

# 타입코드를 정의하지 않으면 기본적으로 정의된 Alphabet()를 사용한다.

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

my_seq = Seq("GATCG", IUPAC.unambiguous_dna)

for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

# 여기에서는 데이터를 객체화하고 그 순번을 출력하는 프로그램이다. for var in enumerate(list)을 사용하면 간략하게 해결된다.

# ================================================================================
from Bio.Seq import Seq

"AAAA".count("AA")
# 2

Seq("AAAA").count("AA")
# 2

또한 각 서열에 포함되어 있는 패턴 서열의 반복수를 구할 수 있다. count()함수를 이용하여 확인해볼 수 있다.

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)
len(my_seq)
# 32

my_seq.count("G")
# 9

# 100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)
# 46.875

# 위의 예제에서는 서열의 GC contents를 구할 수 있다. 

# 이는 전통적인 형태로써 각각 G, C의 개수를 구하여 전체 서열에서 비율을 구하면 된다. 

# ================================================================================
# 이는 다양한 함수 및 수식을 정리하여 사용해야 하지만 아래와 같이 from Bio.SeqUtils import GC를 사용하면 복잡한 수식을 사용하지 않아도 간편하게 값을 확인할 수 있다.

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)

GC(my_seq)
# 46.875

# ================================================================================
# @@ 아래의 예제는 문자열의 슬라이싱 기능을 이용하여 서열을 자르거나 추출하는 것을 보여준다. 

# @@ 서열[x:y:z]처럼 서열을 잘라낼 수 있다. 

# @@ x는 시작위치, y는 마지막 위치, z는 범위를 나타낸다. 

# 아래의 예는 4번 bp 위치에서부터 12번 bp까지 길이 8bp의 서열을 추출하는 것이다.

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)

my_seq[4:12]
# Seq('GATGGGCC', IUPACUnambiguousDNA())

my_seq[0::3]
# Seq('GCTGTAGTAAG', IUPACUnambiguousDNA())

my_seq[1::3]
# Seq('AGGCATGCATC', IUPACUnambiguousDNA())

my_seq[2::3]
# Seq('TAGCTAAGAC', IUPACUnambiguousDNA())

# ================================================================================
# 아래 예제에서는 서열의 처음부터 끝까지 뒤집어서 출력해주는 예제이다.

my_seq[::-1]
# Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG', IUPACUnambiguousDNA())

# ================================================================================
# 또한 서열의 대소문자를 구별하기 위해서는 서열 인스턴스에 .upper()를 사용하면 소문자에서 대문자로 변경할 수 있다.

"GTAC" in dna_seq
# False

"GTAC" in dna_seq.upper()
# True

# ================================================================================
# @@ 또한 서열의 상보서열을 구할 수도 있다. 

# 이 경우에는 .complement(), .reverse_complement()의 기능을 이용하여 상보서열을 아주 쉽게 구할 수 있다.

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
my_seq
# Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPACUnambiguousDNA())

my_seq.complement()
# Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG', IUPACUnambiguousDNA())

my_seq.reverse_complement()
# Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC', IUPACUnambiguousDNA())

# ================================================================================
# 전사 혹은 번역된 서열을 얻기 위해서는 .transcribe(), .translate()를 이용할 수 있다. 

# 전사의 경우에는 T --> U로 DNA에서 mRNA로의 서열로 바꾸어주고 번역의 경우 펩타이트 서열을 보여준다.

coding_dna
# Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

messenger_rna = coding_dna.transcribe()
messenger_rna
# Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

# ================================================================================
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
coding_dna
# Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

coding_dna.translate()
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

# ================================================================================
# @@ 이처럼 Bio.Seq를 이용하여 서열인스턴스를 만들어 간략하게 서열의 변환, 추출, 대소문자, 반복, GC 콘텐츠 등을 구해 볼 수 있다.
