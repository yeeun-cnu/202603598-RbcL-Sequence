
import os

FASTA_PATH=os.path.join("data","processed","rbcl_cleaned.fasta")

def read_fasta(file_path):
    records=[]
    header=""
    sequence=""

    with open(file_path,"r",encoding="utf-8") as file:
        for line in file:
            line=line.strip()

            if line.startswith(">"):
                if header:
                    records.append([header,sequence])

                header=line[1:]
                sequence=""
            else:
                sequence+=line

        if header:
            records.append([header,sequence])

    return records

def main():
    records=read_fasta(FASTA_PATH)

    lengths=[]
    sequence_groups={}

    for header,sequence in records:
        lengths.append(len(sequence))

        if sequence not in sequence_groups:
            sequence_groups[sequence]=[]

        sequence_groups[sequence].append(header)

    print("===== 분석 후보 서열 검사 =====")
    print("전체 서열 수:",len(records))
    print("최단 서열 길이:",min(lengths))
    print("최장 서열 길이:",max(lengths))
    print("평균 서열 길이:",round(sum(lengths)/len(lengths),1))
    print("400 aa 미만:",sum(1 for length in lengths if length<400))


    import csv

    metadata_path=os.path.join("data","metadata","rbcl_landplants_metadata.csv")

    metadata={}

    with open(metadata_path,"r",encoding="utf-8-sig",newline="") as file:
        reader=csv.DictReader(file)

        for row in reader:
            metadata[row["accession"]]=row

    print("\n===== 400 aa 미만 서열 상세 정보 =====")

    for header,sequence in records:
        if len(sequence)<400:
            accession=header.split()[0]
            row=metadata[accession]

            print("Accession:",accession)
            print("생물명:",row["organism"])
            print("길이:",len(sequence))
            print("설명:",row["description"])
            print()


    duplicate_groups=0

    for sequence,headers in sequence_groups.items():
        if len(headers)>1:
            duplicate_groups+=1

    print("동일 서열 그룹 수:",duplicate_groups)
    print("고유 서열 수:",len(sequence_groups))

if __name__=="__main__":
    main()
