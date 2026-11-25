"""
수집한 RbcL FASTA 파일의 기본 품질을 점검하는 프로그램

원본 파일은 수정하거나 삭제하지 않습니다.
서열 길이, 부분 서열 표시, 문자 구성, 완전히 동일한 서열을 확인합니다.
"""

import os


INPUT_PATH = os.path.join("data", "raw", "rbcl_raw.fasta")
VALID_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWYBXZJUO"

# RbcL 전장 후보와 부분 서열을 구분하기 위한 임시 기준입니다.
# 자동 제외 기준이 아니며, 실제 레코드 정보를 추가 확인해야 합니다.
MIN_FULL_LENGTH = 400


def read_fasta(file_path):
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
                records.append({
                    "header": current_header,
                    "sequence": current_sequence
                })

            current_header = line[1:]
            current_sequence = ""

        else:
            current_sequence = current_sequence + line.upper()

    if len(current_header) > 0:
        records.append({
            "header": current_header,
            "sequence": current_sequence
        })

    file.close()
    return records


def check_sequence(sequence):
    invalid_characters = []

    for amino_acid in sequence:
        if amino_acid not in VALID_AMINO_ACIDS:
            if amino_acid not in invalid_characters:
                invalid_characters.append(amino_acid)

    return {
        "length": len(sequence),
        "valid": len(invalid_characters) == 0,
        "invalid_characters": invalid_characters
    }


def find_duplicates(records):
    """완전히 동일한 서열을 묶고 해당 레코드 제목을 반환합니다."""

    sequence_groups = {}

    for record in records:
        sequence = record["sequence"]

        if sequence in sequence_groups:
            sequence_groups[sequence].append(record["header"])
        else:
            sequence_groups[sequence] = [record["header"]]

    duplicate_groups = []

    for sequence in sequence_groups:
        if len(sequence_groups[sequence]) > 1:
            duplicate_groups.append(sequence_groups[sequence])

    return duplicate_groups


def main():
    if os.path.exists(INPUT_PATH) == False:
        print("FASTA 파일을 찾을 수 없습니다.")
        print("먼저 data_collection.py를 실행했는지 확인하세요.")
        return

    records = read_fasta(INPUT_PATH)

    print("===== RbcL FASTA 상세 검사 =====")
    print("전체 서열 수:", len(records))

    if len(records) == 0:
        print("검사할 서열이 없습니다.")
        return

    lengths = []
    invalid_count = 0
    partial_count = 0
    short_count = 0

    for index in range(len(records)):
        record = records[index]
        result = check_sequence(record["sequence"])
        header = record["header"]

        lengths.append(result["length"])

        if result["valid"] == False:
            invalid_count = invalid_count + 1

        is_partial = "partial" in header.lower()
        is_short = result["length"] < MIN_FULL_LENGTH

        if is_partial:
            partial_count = partial_count + 1

        if is_short:
            short_count = short_count + 1

        print()
        print("서열 번호:", index + 1)
        print("제목:", header)
        print("서열 길이:", result["length"])
        print("문자 검사:", result["valid"])
        print("레코드의 partial 표시:", is_partial)
        print("임시 길이 기준 미달:", is_short)

        if result["valid"] == False:
            print("허용되지 않은 문자:", result["invalid_characters"])

    duplicate_groups = find_duplicates(records)

    print()
    print("===== 전체 검사 결과 =====")
    print("전체 서열 수:", len(records))
    print("가장 짧은 서열 길이:", min(lengths))
    print("가장 긴 서열 길이:", max(lengths))
    print("평균 서열 길이:", round(sum(lengths) / len(lengths), 2))
    print("허용되지 않은 문자가 있는 서열 수:", invalid_count)
    print("partial 표시가 있는 서열 수:", partial_count)
    print("길이", MIN_FULL_LENGTH, "미만인 서열 수:", short_count)
    print("완전히 동일한 서열 그룹 수:", len(duplicate_groups))

    if len(duplicate_groups) > 0:
        print()
        print("===== 중복 서열 상세 정보 =====")

        for index in range(len(duplicate_groups)):
            print("중복 그룹", index + 1)

            for header in duplicate_groups[index]:
                print("-", header)

    print()
    print("주의: 길이 기준은 임시 참고 기준이며 자동 제외하지 않습니다.")
    print("원본 FASTA 파일은 변경하지 않았습니다.")


if __name__ == "__main__":
    main()
