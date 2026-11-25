
def needleman_wunsch(seq1,seq2):
    match_score=2
    mismatch_score=-1
    gap_score=-2

    rows=len(seq1)+1
    cols=len(seq2)+1

    matrix=[]

    for i in range(rows):
        matrix.append([0]*cols)

    # 첫 번째 행과 열 초기화
    for i in range(rows):
        matrix[i][0]=i*gap_score

    for j in range(cols):
        matrix[0][j]=j*gap_score

    # 점수 행렬 계산
    for i in range(1,rows):
        for j in range(1,cols):

            if seq1[i-1]==seq2[j-1]:
                diagonal=matrix[i-1][j-1]+match_score
            else:
                diagonal=matrix[i-1][j-1]+mismatch_score

            up=matrix[i-1][j]+gap_score
            left=matrix[i][j-1]+gap_score

            matrix[i][j]=max(diagonal,up,left)

    # 역추적
    aligned1=""
    aligned2=""

    i=len(seq1)
    j=len(seq2)

    while i>0 or j>0:

        if i>0 and j>0:
            if seq1[i-1]==seq2[j-1]:
                score=match_score
            else:
                score=mismatch_score

            if matrix[i][j]==matrix[i-1][j-1]+score:
                aligned1=seq1[i-1]+aligned1
                aligned2=seq2[j-1]+aligned2
                i-=1
                j-=1
                continue

        if i>0 and matrix[i][j]==matrix[i-1][j]+gap_score:
            aligned1=seq1[i-1]+aligned1
            aligned2="-"+aligned2
            i-=1

        else:
            aligned1="-"+aligned1
            aligned2=seq2[j-1]+aligned2
            j-=1

    return aligned1,aligned2,matrix[len(seq1)][len(seq2)]


if __name__=="__main__":

    sequence1="ATCG"
    sequence2="ATG"

    result1,result2,score=needleman_wunsch(sequence1,sequence2)

    print("===== Needleman-Wunsch 정렬 결과 =====")
    print("서열 1:",result1)
    print("서열 2:",result2)
    print("정렬 점수:",score)
