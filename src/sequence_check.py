"""
수집한 RbcL FASTA 파일의 기본 품질을 점검하는 프로그램

이 프로그램은 원본 파일을 수정하거나 서열을 삭제하지 않습니다.
서열 수, 길이, 문자 구성, 중복 여부, 제목에 기록된 생물 정보를 확인합니다.
"""

import os


INPUT_PATH = os.path.join("data", "raw", "rbcl_raw.fasta")

# 표준 아미노산 20종과 일부 단백질 서열에서 사용되는 모호성 문자
VALID_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWYBXZJUO"


def read_fasta(file_path):
    """FASTA 파일을 읽어 제목과 서열을 딕셔너리 목록으로 반환합니다."""

    records = []
    current_header = ""
    current_sequence = ""

    file = open(file_path, "r", encoding="utf-8")

    for line in file:
        line = line.strip()

        if len(line) == 0:
            continue

        if line.startswith(">"):
            if len(current_header) > 0:
                record = {
                    "header": current_header,
                    "sequence": current_sequence
                }
                records.append(record)

            current_header = line[1:]
            current_sequence = ""

        else:
            current_sequence = current_sequence + line.upper()

    if len(current_header) > 0:
        record = {
            "header": current_header,
            "sequence": current_sequence
        }
        records.append(record)

    file.close()

    return records


def check_sequence(sequence):
    """서열의 길이와 아미노산 문자 구성을 검사합니다."""

    invalid_characters = []

    for amino_acid in sequence:
        if amino_acid not in VALID_AMINO_ACIDS:
            if amino_acid not in invalid_characters:
                invalid_characters.append(amino_acid)

    result = {
        "length": len(sequence),
        "valid": len(invalid_characters) == 0,
        "invalid_characters": invalid_characters
    }

    return result


def find_duplicates(records):
    """완전히 동일한 아미노산 서열이 있는지 확인합니다."""

    sequence_counts = {}

    for record in records:
        sequence = record["sequence"]

        if sequence in sequence_counts:
            sequence_counts[sequence] = sequence_counts[sequence] + 1
        else:
            sequence_counts[sequence] = 1

    duplicate_count = 0

    for sequence in sequence_counts:
        if sequence_counts[sequence] > 1:
            duplicate_count = duplicate_count + 1

    return duplicate_count


def main():
    if os.path.exists(INPUT_PATH) == False:
        print("FASTA 파일을 찾을 수 없습니다.")
        print("먼저 data_collection.py를 실행했는지 확인하세요.")
        return

    records = read_fasta(INPUT_PATH)

    print("===== RbcL FASTA 기본 검사 =====")
    print("전체 서열 수:", len(records))

    if len(records) == 0:
        print("검사할 서열이 없습니다.")
        return

    lengths = []
    invalid_count = 0

    for index in range(len(records)):
        record = records[index]
        result = check_sequence(record["sequence"])

        lengths.append(result["length"])

        if result["valid"] == False:
            invalid_count = invalid_count + 1

        print()
        print("서열 번호:", index + 1)
        print("제목:", record["header"])
        print("서열 길이:", result["length"])
        print("문자 검사:", result["valid"])

        if result["valid"] == False:
            print("허용되지 않은 문자:", result["invalid_characters"])

    duplicate_count = find_duplicates(records)

    print()
    print("===== 전체 검사 결과 =====")
    print("전체 서열 수:", len(records))
    print("가장 짧은 서열 길이:", min(lengths))
    print("가장 긴 서열 길이:", max(lengths))
    print("평균 서열 길이:", round(sum(lengths) / len(lengths), 2))
    print("허용되지 않은 문자가 있는 서열 수:", invalid_count)
    print("중복 서열 종류 수:", duplicate_count)
    print()
    print("주의: 이 결과는 기본 검사이며, RbcL 단백질 여부와 생물 종 정보는 추가 검증이 필요합니다.")
    print("원본 FASTA 파일은 변경하지 않았습니다.")


if __name__ == "__main__":
    main()
