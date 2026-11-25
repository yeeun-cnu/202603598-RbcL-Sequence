"""
NCBI에서 식물 RbcL 단백질 서열을 시험 수집하는 프로그램

사용 방법:
1. 프로젝트 최상위 폴더에서 실행합니다.
2. python src/data_collection.py
3. 수집한 FASTA 파일은 data/raw/rbcl_raw.fasta에 저장됩니다.

현재는 연결과 검색 결과를 확인하기 위해 최대 10개를 수집합니다.
시험 수집이 정상적으로 끝난 뒤 수집 규모를 늘립니다.
"""

import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import os
import time


# NCBI Protein 검색 조건
# txid33090은 녹색식물(Viridiplantae) 분류군입니다.
SEARCH_TERM = '(ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit[Protein Name] OR rbcL[Gene Name]) AND txid33090[Organism:exp]'

# 첫 실행에서는 소량만 수집하여 검색 결과를 확인합니다.
MAX_RECORDS = 10

# NCBI 요청에 사용할 기본 주소
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def search_protein_ids():
    """NCBI Protein에서 검색한 레코드의 ID 목록을 반환합니다."""

    parameters = {
        "db": "protein",
        "term": SEARCH_TERM,
        "retmax": MAX_RECORDS,
        "retmode": "xml"
    }

    query_string = urllib.parse.urlencode(parameters)
    url = ESEARCH_URL + "?" + query_string

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "RbcL-Sequence-Analysis/1.0"}
    )

    response = urllib.request.urlopen(request)
    xml_data = response.read()
    response.close()

    root = ET.fromstring(xml_data)

    id_list = []

    for id_element in root.findall(".//Id"):
        id_list.append(id_element.text)

    return id_list


def fetch_fasta(id_list):
    """ID 목록에 해당하는 단백질 서열을 FASTA 형식으로 가져옵니다."""

    if len(id_list) == 0:
        return ""

    parameters = {
        "db": "protein",
        "id": ",".join(id_list),
        "rettype": "fasta",
        "retmode": "text"
    }

    query_string = urllib.parse.urlencode(parameters)
    url = EFETCH_URL + "?" + query_string

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "RbcL-Sequence-Analysis/1.0"}
    )

    response = urllib.request.urlopen(request)
    fasta_text = response.read().decode("utf-8")
    response.close()

    return fasta_text


def save_fasta(fasta_text):
    """가져온 FASTA 데이터를 파일로 저장합니다."""

    output_folder = os.path.join("data", "raw")
    output_path = os.path.join(output_folder, "rbcl_raw.fasta")

    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)

    file = open(output_path, "w", encoding="utf-8")
    file.write(fasta_text)
    file.close()

    return output_path


def main():
    print("NCBI Protein 검색을 시작합니다.")

    try:
        id_list = search_protein_ids()

        print("검색된 레코드 수:", len(id_list))

        if len(id_list) == 0:
            print("검색 결과가 없습니다. 검색 조건을 확인해야 합니다.")
            return

        time.sleep(0.4)

        fasta_text = fetch_fasta(id_list)

        if len(fasta_text) == 0:
            print("FASTA 서열을 가져오지 못했습니다.")
            return

        output_path = save_fasta(fasta_text)

        sequence_count = fasta_text.count(">")

        print("FASTA 서열 수:", sequence_count)
        print("저장 위치:", output_path)
        print("시험 수집이 완료되었습니다.")

    except Exception as error:
        print("데이터 수집 중 오류가 발생했습니다.")
        print("오류 내용:", error)


if __name__ == "__main__":
    main()
