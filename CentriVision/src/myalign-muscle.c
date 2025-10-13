#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char *aligned_seq1;
    char *aligned_seq2;
    int score;
} AlignResult;

AlignResult nw_align(const char *seq1, const char *seq2, int match, int mismatch, int gap) {
    int n = strlen(seq1);
    int m = strlen(seq2);
    int **dp = malloc((n+1) * sizeof(int*));
    for(int i=0;i<=n;i++) dp[i] = malloc((m+1) * sizeof(int));

    // 初始化
    for(int i=0;i<=n;i++) dp[i][0] = i * gap;
    for(int j=0;j<=m;j++) dp[0][j] = j * gap;

    // DP 填表
    for(int i=1;i<=n;i++){
        for(int j=1;j<=m;j++){
            int score_diag = dp[i-1][j-1] + (seq1[i-1]==seq2[j-1]?match:mismatch);
            int score_up = dp[i-1][j] + gap;
            int score_left = dp[i][j-1] + gap;
            int max = score_diag;
            if(score_up>max) max = score_up;
            if(score_left>max) max = score_left;
            dp[i][j] = max;
        }
    }

    // 回溯
    int i=n, j=m, len=0;
    char *a1 = malloc((n+m+1)*sizeof(char));
    char *a2 = malloc((n+m+1)*sizeof(char));
    a1[n+m] = '\0';
    a2[n+m] = '\0';
    int pos = n+m-1;

    while(i>0 || j>0){
        if(i>0 && j>0 && dp[i][j]==dp[i-1][j-1] + (seq1[i-1]==seq2[j-1]?match:mismatch)){
            a1[pos] = seq1[i-1];
            a2[pos] = seq2[j-1];
            i--; j--; pos--;
        } else if(i>0 && dp[i][j]==dp[i-1][j] + gap){
            a1[pos] = seq1[i-1];
            a2[pos] = '-';
            i--; pos--;
        } else {
            a1[pos] = '-';
            a2[pos] = seq2[j-1];
            j--; pos--;
        }
    }

    // 返回结果
    AlignResult res;
    res.aligned_seq1 = strdup(a1+pos+1);
    res.aligned_seq2 = strdup(a2+pos+1);
    res.score = dp[n][m];

    for(int i=0;i<=n;i++) free(dp[i]);
    free(dp); free(a1); free(a2);
    return res;
}

// Python ctypes 调用接口
#include <stdint.h>


// gcc -shared -o myalign.dll myalign.c -Wl,--out-implib,myalign.a
// __declspec(dllexport) AlignResult* nw_align_wrapper(const char* seq1, const char* seq2, int match, int mismatch, int gap){
//     AlignResult* res = malloc(sizeof(AlignResult));
//     *res = nw_align(seq1, seq2, match, mismatch, gap);
//     return res;
// }

// __declspec(dllexport) void free_align_result(AlignResult* res){
//     free(res->aligned_seq1);
//     free(res->aligned_seq2);
//     free(res);
// }


// ==========================
// Linux 修改点：导出符号方式
// Windows: __declspec(dllexport)
// Linux: __attribute__((visibility("default")))
// gcc -fPIC -shared -o myalign.so myalign.c
// ==========================
__attribute__((visibility("default"))) AlignResult* nw_align_wrapper(const char* seq1, const char* seq2, int match, int mismatch, int gap){
    AlignResult* res = malloc(sizeof(AlignResult));
    *res = nw_align(seq1, seq2, match, mismatch, gap);
    return res;
}

__attribute__((visibility("default"))) void free_align_result(AlignResult* res){
    free(res->aligned_seq1);
    free(res->aligned_seq2);
    free(res);
}