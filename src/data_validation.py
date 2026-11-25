"""
FASTA 파일과 메타데이터 CSV가 서로 일치하는지 확인하는 프로그램

사용 방법:
프로젝트 최상위 폴더에서 다음 명령어를 실행합니다.
python src/data_validation.py
"""

import csv
import os


FASTA_PATH = os.path.join("data", "raw", "rbcl_raw.fasta")
CSV_PATH = os.path.join("data", "metadata", "rbcl_metadata.csv")


def read_fasta(path):
    """FASTA 파일에서 accession과 서열 길이를 읽습니다."""

    records = []
    header = ""
    sequence_parts = []

    file = open(path, "r", encoding="utf-8")

    for line in file:
        line = line.strip()

        if len(line) == 0:
            continue

        if line.startswith(">"):
            if len(header) > 0:
                accession = header.split()[0]
                sequence = "".join(sequence_parts)

                records.append({
                    "accession": accession,
                    "description": header,
                    "sequence_length": len(sequence)
                })

            header = line[1:]
            sequence_parts = []

        else:
            sequence_parts.append(line)

    if len(header) > 0:
        accession = header.split()[0]
        sequence = "".join(sequence_parts)

        records.append({
            "accession": accession,
            "description": header,
            "sequence_length": len(sequence)
        })

    file.close()

    return records


def read_metadata(path):
    """CSV 파일에서 메타데이터를 읽습니다."""

    rows = []

    file = open(path, "r", newline="", encoding="utf-8-sig")
    reader = csv.DictReader(file)

    required_columns = [
        "accession",
        "organism",
        "description",
        "sequence_length",
        "partial",
        "taxonomy"
    ]

    if reader.fieldnames is None:
        file.close()
        raise ValueError("CSV 헤더를 읽을 수 없습니다.")

    for column in required_columns:
        if column not in reader.fieldnames:
            file.close()
            raise ValueError("CSV에 필요한 열이 없습니다: " + column)

    for row in reader:
        rows.append(row)

    file.close()

    return rows


def main():
    print("===== FASTA와 CSV 일치 여부 검사 =====")

    if os.path.exists(FASTA_PATH) == False:
        print("FASTA 파일을 찾을 수 없습니다:", FASTA_PATH)
        return

    if os.path.exists(CSV_PATH) == False:
        print("CSV 파일을 찾을 수 없습니다:", CSV_PATH)
        return

    try:
        fasta_records = read_fasta(FASTA_PATH)
        metadata_rows = read_metadata(CSV_PATH)

        print("FASTA 레코드 수:", len(fasta_records))
        print("CSV 데이터 행 수:", len(metadata_rows))

        fasta_accessions = {}
        metadata_accessions = {}

        for record in fasta_records:
            accession = record["accession"]

            if accession in fasta_accessions:
                print("FASTA accession 중복:", accession)

            fasta_accessions[accession] = record

        for row in metadata_rows:
            accession = row["accession"]

            if accession in metadata_accessions:
                print("CSV accession 중복:", accession)

            metadata_accessions[accession] = row

        missing_in_csv = []
        missing_in_fasta = []
        length_mismatches = []
        partial_mismatches = []

        for accession in fasta_accessions:
            if accession not in metadata_accessions:
                missing_in_csv.append(accession)
            else:
                fasta_record = fasta_accessions[accession]
                metadata_row = metadata_accessions[accession]

                try:
                    csv_length = int(metadata_row["sequence_length"])
                except ValueError:
                    csv_length = -1

                if fasta_record["sequence_length"] != csv_length:
                    length_mismatches.append(accession)

                fasta_partial = "partial" in fasta_record["description"].lower()
                csv_partial = metadata_row["partial"].strip().lower() == "true"

                if fasta_partial != csv_partial:
                    partial_mismatches.append(accession)

        for accession in metadata_accessions:
            if accession not in fasta_accessions:
                missing_in_fasta.append(accession)

        print("CSV에 없는 FASTA accession 수:", len(missing_in_csv))
        print("FASTA에 없는 CSV accession 수:", len(missing_in_fasta))
        print("서열 길이 불일치 수:", len(length_mismatches))
        print("partial 표시 불일치 수:", len(partial_mismatches))

        if len(missing_in_csv) == 0 and len(missing_in_fasta) == 0 and len(length_mismatches) == 0 and len(partial_mismatches) == 0 and len(fasta_records) == len(metadata_rows):
            print("검사 결과: FASTA와 CSV가 서로 일치합니다.")
        else:
            print("검사 결과: 불일치 항목을 확인해야 합니다.")

        print("원본 파일은 변경하지 않았습니다.")

    except Exception as error:
        print("검사 중 오류가 발생했습니다.")
        print("오류 내용:", error)


if __name__ == "__main__":
    main()
