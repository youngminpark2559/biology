
# Biopython에서 서열 레코드 객체는 서열 그 자체에 대한 정보 뿐만 아니라 
# 서열의 이름, description, 일부의 경우 annotation, sub feature 정보를 포함한다 (Biopython/Bio.SeqRecord 참고).

# SeqFeature 객체를 이용하여 서열 레코드 객체에 gene, CDS 등 새로운 feature 정보를 추가할 수 있는데 다음과 같은 경우 유용하게 이용될 수 있다.

# GFF 파일과 FASTA 파일을 이용하여 Genbank 파일을 생성하는 경우
# augustus는 대표적인 gene prediction tool로써 유전자 구조 예측을 위하여 종 특이적인 매트릭스를 필요로 한다.

# 대표적인 종의 경우 기구축된 매트릭스를 제공하고 있으나 de novo 종의 경우 사용자가 직접 생성할 수 있도록 etraining이라는 프로그램을 별도로 제공한다.

# 이때 etraining의 입력 파일이 Genbank 파일 포맷이어야 한다.

# RNAseq이나 Reference Protein을 genome에 맵핑 후 full length로 얻어지는 evidence gene model (GFF) 정보를 이용하여 training을 위한 Genbank 파일을 만들 수 있다.

# 이미 만들어진 Genbank 파일에 새로운 feature (gene, CDS, mRNA, UTR 등) 정보를 추가하는 경우 SeqFeature 객체를 이용하기 위해서는 먼저 다음과 같이 객체를 불러와야 한다.

# FeatureLocation 객체는 각 feature의 위치 정보를 입력할 때 필요하므로 함께 호출되도록 한다.

from Bio.SeqFeature import SeqFeature, FeatureLocation

# SeqFeature는 서열 격체에 정보를 추가하는 것이므로 서열 객체의 존재가 필수적이다. 

# "genbank"라고 이름지어진 서열 객체에 feature 정보를 입력하는 방법은 다음과 같다. 

# 참고로 augustus training을 위한 Genbank 파일 생성시에는 source feature를 먼저 추가해주어야 한다.

feature = SeqFeature(FeatureLocation(start=0, end=50)), type='source')

# 생성된 feature를 genbank 서열 객체에 추가한다.
genbank.features.append(feature)

# gene feature를 생성하고 genbank 서열 객체에 추가해보자. 이때 각 feature의 방향성 (strand) 정보가 같이 입력되어야 하며, 정방향일 경우 +1, reverse complement일 경우 -1을 입력해야 한다.
SeqFeature(FeatureLocation(start=5, end=30), type='gene', strand=-1)
genbank.features.append(feature)

# 하나의 유전자에 여러 CDS가 존재하는 경우 각각의 FeatureLocation을 생성 후 병합하여 하나의 SeqFeature로 추가할 수 있다.

for index, (sp,ep) in enumerate(cds_pos):
    if index == 0:
        combined = FeatureLocation(sp-1, ep)
    else:
        combined = combined + FeatureLocation(sp-1, ep)

feature = SeqFeature(combined, type='CDS', strand=-1, qualifiers={'translate':'MVNVPKERRTYC'})
genbank.features.append(feature)

# ================================================================================
# 이렇게 만들어진 Genbank 서열 객체는 Biopython/Bio.SeqIO 객체를 이용하여 파일에 쓸 수 있다 (output format은 'genbank'로 설정한다).

SeqIO.write(genbank, output, 'genbank')

# Genbank 파일은 비교적 간단한 FASTA 파일과 달리 다양한 정보를 포함하고 있고 들여쓰기 및 공백의 유무가 중요하게 작용한다. 

# Python과 같은 프로그래밍 언어에 능통하다면 굳이 Biopython을 이용하지 않고도 충분히 생성이 가능하지만 Biopython을 이용한다면 실행 오류와 코딩 라인을 줄일 수 있어 훨씬 유용할 것이다.