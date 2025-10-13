// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>

// typedef struct {
//     char* aligned_seq1;
//     char* aligned_seq2;
//     int score;
// } AlignResult;

// int max3(int a,int b,int c){ return (a>b?(a>c?a:c):(b>c?b:c)); }

// #define IDX(i,j,len2) ((i)*(len2+1)+(j))

// AlignResult* nw_align_iterative(const char* seq1,const char* seq2,int match,int mismatch,int gap_open,int gap_extend){
//     int len1 = strlen(seq1);
//     int len2 = strlen(seq2);

//     int* dp = (int*)malloc((len1+1)*(len2+1)*sizeof(int));
//     char* traceback = (char*)malloc((len1+1)*(len2+1)*sizeof(char));

//     dp[IDX(0,0,len2)] = 0;
//     traceback[IDX(0,0,len2)] = '0';

//     for(int i=1;i<=len1;i++){
//         dp[IDX(i,0,len2)] = gap_open + (i-1)*gap_extend;
//         traceback[IDX(i,0,len2)] = 'U';
//     }
//     for(int j=1;j<=len2;j++){
//         dp[IDX(0,j,len2)] = gap_open + (j-1)*gap_extend;
//         traceback[IDX(0,j,len2)] = 'L';
//     }

//     for(int i=1;i<=len1;i++){
//         for(int j=1;j<=len2;j++){
//             int score_diag = dp[IDX(i-1,j-1,len2)] + (seq1[i-1]==seq2[j-1]?match:mismatch);
//             int score_up   = dp[IDX(i-1,j,len2)] + ((traceback[IDX(i-1,j,len2)]=='U')?gap_extend:gap_open);
//             int score_left = dp[IDX(i,j-1,len2)] + ((traceback[IDX(i,j-1,len2)]=='L')?gap_extend:gap_open);
//             dp[IDX(i,j,len2)] = max3(score_diag,score_up,score_left);

//             if(dp[IDX(i,j,len2)]==score_diag) traceback[IDX(i,j,len2)]='D';
//             else if(dp[IDX(i,j,len2)]==score_up) traceback[IDX(i,j,len2)]='U';
//             else traceback[IDX(i,j,len2)]='L';
//         }
//     }

//     // 回溯
//     int i=len1,j=len2;
//     int alloc_size = len1+len2+2;
//     char* aln1 = (char*)malloc(alloc_size);
//     char* aln2 = (char*)malloc(alloc_size);
//     aln1[alloc_size-1]='\0';
//     aln2[alloc_size-1]='\0';
//     int p = alloc_size-2;

//     while(i>0 || j>0){
//         if(i>0 && j>0 && traceback[IDX(i,j,len2)]=='D'){
//             aln1[p]=seq1[i-1]; aln2[p]=seq2[j-1]; i--; j--; p--;
//         } else if(i>0 && traceback[IDX(i,j,len2)]=='U'){
//             aln1[p]=seq1[i-1]; aln2[p]='-'; i--; p--;
//         } else { // L
//             aln1[p]='-'; aln2[p]=seq2[j-1]; j--; p--;
//         }
//     }

//     char* final1 = strdup(aln1+p+1);
//     char* final2 = strdup(aln2+p+1);

//     // 计算分数
//     int score=0;
//     for(int k=0;k<strlen(final1);k++){
//         if(final1[k]==final2[k]) score+=match;
//         else if(final1[k]=='-'||final2[k]=='-') score+=gap_extend;
//         else score+=mismatch;
//     }

//     free(aln1); free(aln2); free(dp); free(traceback);

//     AlignResult* res = (AlignResult*)malloc(sizeof(AlignResult));
//     res->aligned_seq1 = final1;
//     res->aligned_seq2 = final2;
//     res->score = score;
//     return res;
// }

// void free_align_result(AlignResult* res){
//     if(res){
//         if(res->aligned_seq1) free(res->aligned_seq1);
//         if(res->aligned_seq2) free(res->aligned_seq2);
//         free(res);
//     }
// }

// AlignResult* nw_align_wrapper(const char* seq1,const char* seq2,int match,int mismatch,int gap){
//     int gap_open=gap;
//     int gap_extend=gap/2;
//     return nw_align_iterative(seq1,seq2,match,mismatch,gap_open,gap_extend);
// }


