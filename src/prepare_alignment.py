
import os

INPUT_PATH=os.path.join("data","processed","rbcl_cleaned.fasta")
OUTPUT_PATH=os.path.join("data","processed","rbcl_alignment_input.fasta")

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
    records=read_fasta(INPUT_PATH)

    selected=[]
    excluded=[]

    for header,sequence in records:
        if len(sequence)>=400:
            selected.append([header,sequence])
        else:
            excluded.append([header,sequence])

    with open(OUTPUT_PATH,"w",encoding="utf-8") as file:
        for header,sequence in selected:
            file.write(">"+header+"\n")
            file.write(sequence+"\n")

    print("===== 최종 정렬 후보 준비 =====")
    print("기존 분석 후보:",len(records))
    print("정렬 후보:",len(selected))
    print("길이 기준 제외:",len(excluded))
    print("저장 위치:",OUTPUT_PATH)

if __name__=="__main__":
    main()
