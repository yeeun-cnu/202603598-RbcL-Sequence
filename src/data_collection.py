"""
NCBI에서 식물 RbcL 단백질 서열과 메타데이터를 수집하는 프로그램

사용 방법:
1. 프로젝트 최상위 폴더에서 실행합니다.
2. python src/data_collection.py
3. FASTA 파일은 data/raw/rbcl_raw.fasta에 저장됩니다.
4. 메타데이터 CSV는 data/metadata/rbcl_metadata.csv에 저장됩니다.

현재는 기능 확인을 위해 최대 10개 레코드만 수집합니다.
"""

import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import csv
import os
import time


# NCBI Protein 검색 조건
# txid33090은 녹색식물(Viridiplantae) 분류군입니다.
SEARCH_TERM = '(ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit[Protein Name] OR rbcL[Gene Name]) AND txid33090[Organism:exp]'

# 첫 실행에서는 소량만 수집하여 결과를 확인합니다.
MAX_RECORDS = 100

# NCBI 요청 주소
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

USER_AGENT = "RbcL-Sequence-Analysis/1.0"


def send_request(url):
    """NCBI에 요청을 보내고 응답 내용을 반환합니다."""

    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT}
    )

    response = urllib.request.urlopen(request, timeout=60)
    data = response.read()
    response.close()

    return data


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

    xml_data = send_request(url)
    root = ET.fromstring(xml_data)

    id_list = []

    for id_element in root.findall(".//Id"):
        if id_element.text:
            id_list.append(id_element.text)

    return id_list


def fetch_protein_records(id_list):
    """ID 목록에 해당하는 단백질 GenBank XML 레코드를 가져옵니다."""

    if len(id_list) == 0:
        return None

    parameters = {
        "db": "protein",
        "id": ",".join(id_list),
        "rettype": "gb",
        "retmode": "xml"
    }

    query_string = urllib.parse.urlencode(parameters)
    url = EFETCH_URL + "?" + query_string

    xml_data = send_request(url)
    root = ET.fromstring(xml_data)

    return root


def get_record_text(record, tag_name):
    """XML 레코드에서 지정한 항목의 텍스트를 가져옵니다."""

    element = record.find(tag_name)

    if element is None:
        return ""

    if element.text is None:
        return ""

    return element.text.strip()


def parse_records(root):
    """XML 레코드를 FASTA 데이터와 메타데이터 목록으로 변환합니다."""

    fasta_records = []
    metadata_rows = []

    for record in root.findall(".//GBSeq"):
        accession = get_record_text(record, "GBSeq_accession-version")
        description = get_record_text(record, "GBSeq_definition")
        organism = get_record_text(record, "GBSeq_organism")
        taxonomy = get_record_text(record, "GBSeq_taxonomy")
        sequence = get_record_text(record, "GBSeq_sequence")

        sequence = sequence.upper()
        sequence_length = len(sequence)

        if len(accession) == 0 or len(sequence) == 0:
            continue

        partial = "partial" in description.lower()

        fasta_header = ">" + accession + " " + description
        fasta_record = fasta_header + "\n" + sequence

        fasta_records.append(fasta_record)

        metadata_row = [
            accession,
            organism,
            description,
            sequence_length,
            partial,
            taxonomy
        ]

        metadata_rows.append(metadata_row)

    return fasta_records, metadata_rows


def save_fasta(fasta_records):
    """FASTA 레코드 목록을 파일로 저장합니다."""

    output_folder = os.path.join("data", "raw")
    output_path = os.path.join(output_folder, "rbcl_raw.fasta")

    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)

    file = open(output_path, "w", encoding="utf-8")
    file.write("\n".join(fasta_records))

    if len(fasta_records) > 0:
        file.write("\n")

    file.close()

    return output_path


def save_metadata(metadata_rows):
    """메타데이터 목록을 CSV 파일로 저장합니다."""

    output_folder = os.path.join("data", "metadata")
    output_path = os.path.join(output_folder, "rbcl_metadata.csv")

    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)

    file = open(output_path, "w", newline="", encoding="utf-8-sig")
    writer = csv.writer(file)

    writer.writerow([
        "accession",
        "organism",
        "description",
        "sequence_length",
        "partial",
        "taxonomy"
    ])

    for row in metadata_rows:
        writer.writerow(row)

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

        # NCBI 요청 사이에 잠시 대기합니다.
        time.sleep(0.4)

        root = fetch_protein_records(id_list)

        if root is None:
            print("가져올 레코드가 없습니다.")
            return

        fasta_records, metadata_rows = parse_records(root)

        if len(fasta_records) == 0:
            print("유효한 단백질 서열을 가져오지 못했습니다.")
            return

        fasta_path = save_fasta(fasta_records)
        metadata_path = save_metadata(metadata_rows)

        print("FASTA 서열 수:", len(fasta_records))
        print("메타데이터 레코드 수:", len(metadata_rows))
        print("FASTA 저장 위치:", fasta_path)
        print("CSV 저장 위치:", metadata_path)
        print("데이터 수집이 완료되었습니다.")

    except Exception as error:
        print("데이터 수집 중 오류가 발생했습니다.")
        print("오류 내용:", error)


if __name__ == "__main__":
    main()
