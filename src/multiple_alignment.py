
import os
from pairwise_alignment import needleman_wunsch

INPUT_PATH=os.path.join("data","processed","rbcl_alignment_input.fasta")
OUTPUT_PATH=os.path.join("data","processed","rbcl_aligned.fasta")

# 최초 검증은 10개, 성공 후 60개로 변경
LIMIT=60


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


def alignment_score(seq1,seq2):
    match_score=2
    mismatch_score=-1
    gap_score=-2

    previous=[j*gap_score for j in range(len(seq2)+1)]

    for i in range(1,len(seq1)+1):
        current=[i*gap_score]

        for j in range(1,len(seq2)+1):
            if seq1[i-1]==seq2[j-1]:
                diagonal=previous[j-1]+match_score
            else:
                diagonal=previous[j-1]+mismatch_score

            up=previous[j]+gap_score
            left=current[j-1]+gap_score

            current.append(max(diagonal,up,left))

        previous=current

    return previous[-1]


def find_center(records):
    sequences=[record[1] for record in records]
    total_scores=[]

    for i in range(len(sequences)):
        total=0

        for j in range(len(sequences)):
            if i!=j:
                total+=alignment_score(sequences[i],sequences[j])

        total_scores.append(total)

    center_index=total_scores.index(max(total_scores))

    return center_index,total_scores



def extract_insertions(center,other):
    center_length=sum(1 for char in center if char!="-")

    insertions=["" for i in range(center_length+1)]
    aligned_other=[]
    position=0

    for a,b in zip(center,other):
        if a=="-":
            insertions[position]+=b
        else:
            aligned_other.append(b)
            position+=1

    return insertions,aligned_other


def center_star(records):
    center_index,scores=find_center(records)

    center_sequence=records[center_index][1]

    pair_data=[]
    max_insertions=[0 for i in range(len(center_sequence)+1)]

    for index,record in enumerate(records):
        if index==center_index:
            continue

        aligned_center,aligned_other,score=needleman_wunsch(
            center_sequence,record[1]
        )

        insertions,residues=extract_insertions(
            aligned_center,aligned_other
        )

        for i in range(len(insertions)):
            if len(insertions[i])>max_insertions[i]:
                max_insertions[i]=len(insertions[i])

        pair_data.append([index,insertions,residues])


    aligned_records=[None for i in range(len(records))]

    center_aligned=""

    for i in range(len(center_sequence)+1):
        center_aligned+="-"*max_insertions[i]

        if i<len(center_sequence):
            center_aligned+=center_sequence[i]

    aligned_records[center_index]=[
        records[center_index][0],
        center_aligned
    ]

    for index,insertions,residues in pair_data:
        aligned_sequence=""

        for i in range(len(center_sequence)+1):
            padding=max_insertions[i]-len(insertions[i])

            aligned_sequence+=insertions[i]
            aligned_sequence+="-"*padding

            if i<len(center_sequence):
                aligned_sequence+=residues[i]

        aligned_records[index]=[
            records[index][0],
            aligned_sequence
        ]

    return aligned_records,center_index,scores



def main():
    records=read_fasta(INPUT_PATH)
    records=records[:LIMIT]

    if len(records)<2:
        print("정렬하려면 최소 2개의 서열이 필요합니다.")
        return

    print("===== Center-Star 다중서열정렬 =====")
    print("입력 서열 수:",len(records))

    aligned_records,center_index,scores=center_star(records)

    print("\n===== 정렬 품질 기초 검사 =====")

    original_lengths=[len(sequence) for header,sequence in records]

    gap_count=0
    gap_sequences=0

    for header,sequence in aligned_records:
        count=sequence.count("-")
        gap_count+=count

        if count>0:
            gap_sequences+=1

    print("정렬 전 최단 길이:",min(original_lengths))
    print("정렬 전 최장 길이:",max(original_lengths))
    print("정렬 후 gap 총 개수:",gap_count)
    print("Gap이 포함된 서열 수:",gap_sequences)


    aligned_records,center_index,scores=center_star(records)

    print("\n===== 원본 서열 보존 검사 =====")

    preserved_count=0

    for i in range(len(records)):
        original=records[i][1]
        aligned=aligned_records[i][1]

        restored=aligned.replace("-","")

        if original==restored:
            preserved_count+=1
        else:
            print("불일치:",records[i][0])

    print("원본과 일치하는 서열:",preserved_count)
    print("원본 서열 보존 여부:",preserved_count==len(records))

    lengths=[]

    for header,sequence in aligned_records:
        lengths.append(len(sequence))


    lengths=[]

    for header,sequence in aligned_records:
        lengths.append(len(sequence))

    with open(OUTPUT_PATH,"w",encoding="utf-8") as file:
        for header,sequence in aligned_records:
            file.write(">"+header+"\n")
            file.write(sequence+"\n")

    print("중심 서열:",records[center_index][0])
    print("정렬된 서열 수:",len(aligned_records))
    print("정렬 후 서열 길이:",min(lengths))
    print("정렬 후 서열 길이:",max(lengths))
    print("모든 서열 길이 동일:",len(set(lengths))==1)
    print("저장 위치:",OUTPUT_PATH)


if __name__=="__main__":
    main()
