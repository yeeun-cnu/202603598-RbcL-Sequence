
import csv
import os

CSV_PATH=os.path.join("data","metadata","rbcl_landplants_metadata.csv")
FASTA_PATH=os.path.join("data","raw","rbcl_landplants_raw.fasta")

OUTPUT_PATH=os.path.join("data","processed","rbcl_cleaned.fasta")
REPORT_PATH=os.path.join("data","processed","rbcl_cleaning_report.csv")

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
    os.makedirs("data/processed",exist_ok=True)

    fasta_records=read_fasta(FASTA_PATH)

    metadata={}

    with open(CSV_PATH,"r",encoding="utf-8-sig",newline="") as file:
        reader=csv.DictReader(file)

        for row in reader:
            metadata[row["accession"]]=row

    selected=[]
    report=[]

    for header,sequence in fasta_records:
        accession=header.split()[0]
        row=metadata[accession]

        if row["partial"].lower()=="true":
            status="Excluded: partial"
        else:
            status="Selected"

            if len(sequence)<400:
                status="Selected: review length"

            selected.append([header,sequence])

        report.append({
            "accession":accession,
            "organism":row["organism"],
            "sequence_length":len(sequence),
            "partial":row["partial"],
            "status":status
        })

    with open(OUTPUT_PATH,"w",encoding="utf-8") as file:
        for header,sequence in selected:
            file.write(">"+header+"\n")
            file.write(sequence+"\n")

    with open(REPORT_PATH,"w",encoding="utf-8-sig",newline="") as file:
        fields=["accession","organism","sequence_length","partial","status"]
        writer=csv.DictWriter(file,fieldnames=fields)
        writer.writeheader()
        writer.writerows(report)

    print("===== 데이터 정제 결과 =====")
    print("전체 서열 수:",len(fasta_records))
    print("분석 후보 수:",len(selected))
    print("제외된 서열 수:",len(fasta_records)-len(selected))
    print("분석용 FASTA:",OUTPUT_PATH)
    print("정제 기록 CSV:",REPORT_PATH)

if __name__=="__main__":
    main()
