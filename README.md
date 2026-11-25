# 202603598-RbcL-Sequence

NCBI Protein 데이터베이스에서 육상식물(Embryophyta)의 RbcL 단백질 서열을 수집하고, 품질 관리와 직접 구현한 서열 정렬 알고리즘을 통해 아미노산 보존성 및 식물 간 서열 차이를 분석하는 프로젝트입니다.

## 1. 프로젝트 개요

RbcL(Ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit)은 광합성의 탄소 고정 과정에 관여하는 루비스코(Rubisco)의 대단위체입니다. 본 프로젝트에서는 특정 소수의 식물 종을 미리 정하는 대신, NCBI Protein 데이터베이스에서 확보할 수 있는 육상식물의 RbcL 단백질 서열을 수집합니다.

수집한 서열의 품질을 검사하고 분석에 적합한 서열을 선별한 뒤, Needleman-Wunsch 쌍별 정렬과 Center-Star 다중 서열 정렬을 직접 구현하여 서열 간 공통점과 차이를 비교합니다. 이후 정렬 결과를 바탕으로 위치별 아미노산 보존성과 분류군별 유사도를 분석하는 것을 목표로 합니다.

## 2. 연구 목적 및 질문

- NCBI Protein에서 육상식물 RbcL 단백질 서열을 얼마나 확보할 수 있는가?
- 수집한 서열 중 불완전하거나 분석에 적합하지 않은 서열은 얼마나 되는가?
- 다양한 식물의 RbcL 단백질 서열에서 공통적으로 보존되는 아미노산 위치는 어디인가?
- 변이가 자주 나타나는 위치는 어디이며, 분류군에 따라 차이가 있는가?
- 직접 구현한 정렬 알고리즘의 결과는 기존 정렬 도구의 결과와 어떤 차이가 있는가?

## 3. 분석 대상 및 데이터 출처

- 분석 대상: 육상식물(Embryophyta)의 RbcL 단백질 서열
- 데이터베이스: NCBI Protein
- 검색 범위: NCBI Taxonomy의 Embryophyta 분류군(taxonomy ID: 3193)을 기준으로 검색
- 초기 수집 규모: 단백질 서열 레코드 100개
- 데이터 형식: FASTA 및 CSV

여기서 100개는 종의 수가 아니라 수집한 단백질 서열 레코드의 수입니다. 동일한 생물종에서 여러 레코드가 수집될 수 있으므로, 레코드 수와 종의 수를 동일하게 해석하지 않습니다.

메타데이터에는 accession, 생물종 정보, 레코드 설명, 서열 길이, partial 여부 및 분류학적 정보 등을 기록합니다.

## 4. 데이터 처리 및 품질 관리

데이터는 원본, 메타데이터, 처리 결과를 구분하여 저장합니다.

1. NCBI Protein에서 RbcL 단백질 서열과 관련 정보를 수집합니다.
2. FASTA와 CSV의 레코드 수, accession, 서열 길이 및 partial 표시가 일치하는지 확인합니다.
3. 아미노산 문자와 서열 길이 등을 검사합니다.
4. partial로 표시된 서열을 정제 단계에서 제외합니다.
5. 정제된 서열 중 길이가 400개 아미노산 미만인 서열은 다중 서열 정렬 대상에서 제외합니다.
6. 원본 파일은 보존하고 정제 결과와 정렬 입력 파일을 별도로 생성합니다.

### 현재 확인된 데이터 처리 결과

| 처리 단계 | 서열 수 | 설명 |
|---|---:|---|
| 최초 수집 | 100개 | 육상식물 RbcL 단백질 레코드 |
| partial 제외 후 | 63개 | partial로 표시된 37개 제외 |
| 길이 기준 적용 후 | 60개 | 400 aa 미만인 3개 추가 제외 |
| 다중 서열 정렬 | 60개 | Center-Star 정렬 입력 및 결과 |

FASTA와 CSV 사이의 accession, 서열 길이 및 partial 표시 일치 여부 검사에서는 불일치가 발견되지 않았습니다. 정제 단계에서는 partial 서열 37개를 제외했으며, 정렬 준비 단계에서는 400 aa 미만인 서열 3개를 추가로 제외했습니다.

## 5. 분석 방법

### 5.1. 데이터 수집 및 검증

`src/data_collection.py`를 이용해 NCBI에서 서열과 메타데이터를 수집합니다. `src/data_validation.py`는 FASTA와 CSV의 레코드 대응 관계 및 주요 정보의 일치 여부를 확인합니다. `src/sequence_check.py`는 원본 서열의 길이, 문자 구성, partial 여부 및 중복 등을 점검합니다.

### 5.2. 서열 정제 및 정렬 준비

`src/data_cleaning.py`는 partial 서열을 제외한 FASTA와 정제 보고서를 생성합니다. `src/cleaned_data_check.py`로 정제된 서열을 다시 검사하고, `src/prepare_alignment.py`에서 길이 기준을 적용하여 다중 서열 정렬 입력 파일을 준비합니다.

### 5.3. 쌍별 서열 정렬

`src/pairwise_alignment.py`에 Needleman-Wunsch 전역 정렬 알고리즘을 구현했습니다. 현재 사용한 점수 설정은 다음과 같습니다.

- 일치(match): +2
- 불일치(mismatch): -1
- Gap penalty: -2

이는 교육 및 알고리즘 구현을 위한 단순 점수 체계이며, 단백질 전용 치환 행렬을 적용한 결과와는 차이가 있을 수 있습니다.

### 5.4. 다중 서열 정렬

`src/multiple_alignment.py`에 Center-Star 방식의 다중 서열 정렬(Multiple Sequence Alignment, MSA)을 구현했습니다. 쌍별 정렬을 이용해 중심 서열을 정하고, 나머지 서열을 중심 서열에 맞춰 정렬합니다.

현재 실행 결과는 다음과 같습니다.

- 정렬 입력 서열: 60개
- 정렬 후 서열 길이: 492 aa columns
- 정렬 후 전체 gap 개수: 947개
- Gap이 포함된 서열: 60개
- Gap 제거 후 원본 서열과 일치한 서열: 60개
- 모든 정렬 서열의 길이 동일 여부: True

원본 서열 보존 검사는 정렬 결과에서 gap을 제거했을 때 입력 서열이 복원되는지 확인한 것입니다. 이 검사는 서열이 손실되거나 순서가 잘못 바뀌지 않았는지 확인하지만, 정렬 위치 자체가 생물학적으로 최적이라는 것을 보장하지는 않습니다.

## 6. 프로젝트 구조

```text
202603598-RbcL-Sequence/
├── data/
│   ├── raw/
│   │   ├── rbcl_landplants_raw.fasta
│   │   └── rbcl_raw.fasta
│   ├── metadata/
│   │   ├── rbcl_landplants_metadata.csv
│   │   └── rbcl_metadata.csv
│   └── processed/
│       ├── rbcl_cleaned.fasta
│       ├── rbcl_cleaning_report.csv
│       ├── rbcl_alignment_input.fasta
│       └── rbcl_aligned.fasta
├── src/
│   ├── data_collection.py
│   ├── sequence_check.py
│   ├── data_validation.py
│   ├── data_cleaning.py
│   ├── cleaned_data_check.py
│   ├── prepare_alignment.py
│   ├── pairwise_alignment.py
│   └── multiple_alignment.py
└── README.md
```

## 7. 실행 순서

프로젝트 최상위 폴더에서 아래 명령을 순서대로 실행합니다.

```powershell
python src/data_collection.py
python src/sequence_check.py
python src/data_validation.py
python src/data_cleaning.py
python src/cleaned_data_check.py
python src/prepare_alignment.py
python src/pairwise_alignment.py
python src/multiple_alignment.py
```

각 스크립트의 입력 및 출력 경로는 코드에 정의되어 있습니다. 데이터 수집 스크립트는 NCBI에 접속하므로 네트워크 연결이 필요합니다. 이미 수집한 데이터를 검사하거나 정렬하는 단계에서는 해당 입력 파일이 프로젝트 폴더에 있어야 합니다.

## 8. 현재 진행 상황

- [x] 연구 주제 및 분석 방향 설정
- [x] 육상식물 RbcL 단백질 서열 100개 수집
- [x] FASTA와 CSV 메타데이터 일치 여부 검증
- [x] 원본 서열의 길이, 문자 구성, partial 및 중복 정보 검사
- [x] partial 서열 37개 제외 및 정제 데이터 생성
- [x] 정제 데이터 검사
- [x] 400 aa 길이 기준으로 정렬 대상 60개 선정
- [x] Needleman-Wunsch 쌍별 정렬 알고리즘 구현
- [x] Center-Star 다중 서열 정렬 알고리즘 구현
- [x] 정렬 결과의 길이 일치 및 원본 서열 보존 여부 검사
- [ ] MAFFT 또는 Clustal Omega 결과와 비교
- [ ] 위치별 아미노산 보존성 및 Shannon entropy 분석
- [ ] 종간·분류군별 유사도 분석
- [ ] 결과 시각화 및 생물학적 해석

## 9. 해석상의 유의 사항 및 향후 검증

Center-Star 방식은 여러 서열을 정렬하기 위한 실용적인 접근법이지만, 가능한 모든 다중 서열 정렬 중 전역적으로 최적의 결과를 보장하지는 않습니다. 또한 현재 구현은 단순 match/mismatch/gap 점수 체계를 사용하므로 단백질 진화적 치환을 충분히 반영하지 못할 수 있습니다.

정렬 결과에서 gap을 제거했을 때 원본 서열이 복원되고 모든 정렬 서열의 길이가 같다는 것은 기본적인 구조 검사가 통과했다는 의미입니다. 정렬의 생물학적 타당성을 평가하려면 MAFFT 또는 Clustal Omega 등 기존 도구의 결과와 비교하고, 보존 영역 및 알려진 RbcL 특징과 대조하는 추가 검증이 필요합니다.

또한 NCBI 데이터베이스의 서열은 특정 분류군이나 연구 대상에 편중될 수 있으므로, 수집한 100개 레코드가 육상식물 전체의 다양성을 균등하게 대표한다고 볼 수 없습니다. 분석 결과는 사용한 데이터의 범위와 품질 관리 기준을 고려하여 해석해야 합니다.

## 10. 참고 데이터베이스 및 도구

- [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/)
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [MAFFT](https://mafft.cbrc.jp/alignment/server/)
- [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/)
