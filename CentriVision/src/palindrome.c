/*#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

char complement(char base) {
    switch (toupper(base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default:  return base;
    }
}

// 检查回文类型，返回值：1=EQ, 2=PS, 3=RC, 0=非回文
int check_palindrome(const char* seq1, const char* seq2, int length) {
    int eq = 1, ps = 1, rc = 1;
    for (int i = 0; i < length; i++) {
        char a = toupper(seq1[i]);
        char b = toupper(seq2[i]);
        if (a != b) eq = 0;
        if (a != toupper(seq2[length - 1 - i])) ps = 0;
        if (a != complement(seq2[i])) rc = 0;
    }
    if (eq) return 1;
    if (ps) return 2;
    if (rc) return 3;
    return 0;
}

// 全序列扫描回文
int find_palindromes(const char* seq, int seqlength, int k, int reach, int* results, int max_results) {
    int count = 0;
    for (int i = 0; i < seqlength - 2*k; i += k) {
        const char* seq1 = seq + i;
        for (int j = i + 1; j < i + reach && j < seqlength - k; j++) {
            const char* seq2 = seq + j;
            int t = check_palindrome(seq1, seq2, k);
            if (t > 0) {
                if (count + 2 > max_results) return count; // 避免越界
                // 第一条
                results[count*5 + 0] = i;
                results[count*5 + 1] = i + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = j;
                count++;
                // 第二条
                results[count*5 + 0] = j;
                results[count*5 + 1] = j + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = i;
                count++;
            }
        }
    }
    return count;
}

// 分块扫描回文
int find_palindromes_chunk(const char* seq, int seqlength, int k, int reach, int start_idx, int end_idx, int* results) {
    int count = 0;
    if (start_idx < 0) start_idx = 0;
    if (end_idx > seqlength) end_idx = seqlength;
    for (int i = start_idx; i < end_idx - k; i += k) {
        const char* seq1 = seq + i;
        for (int j = i + 1; j < i + reach && j < end_idx - k; j++) {
            const char* seq2 = seq + j;
            int t = check_palindrome(seq1, seq2, k);
            if (t > 0) {
                // 直接写入数组
                results[count*5 + 0] = i;
                results[count*5 + 1] = i + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = j;
                count++;
                results[count*5 + 0] = j;
                results[count*5 + 1] = j + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = i;
                count++;
            }
        }
    }
    return count;
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

char complement(char base) {
    switch (toupper(base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default:  return base;
    }
}

// 检查序列是否包含'N'或所有碱基相同
int contains_N_or_homogeneous(const char* seq, int length) {
    char first = toupper(seq[0]);
    for (int i = 0; i < length; i++) {
        char base = toupper(seq[i]);
        if (base == 'N') return 1;
        if (base != first) return 0;
    }
    return 1;
}

// 检查回文类型，返回值：1=EQ, 2=PS, 3=RC, 0=非回文
int check_palindrome(const char* seq1, const char* seq2, int length) {
    int eq = 1, ps = 1, rc = 1;
    for (int i = 0; i < length; i++) {
        char a = toupper(seq1[i]);
        char b = toupper(seq2[i]);
        if (a != b) eq = 0;
        if (a != toupper(seq2[length - 1 - i])) ps = 0;
        if (a != complement(seq2[i])) rc = 0;
    }
    if (eq) return 1;
    if (ps) return 2;
    if (rc) return 3;
    return 0;
}

// 全序列扫描回文（可能未使用）
int find_palindromes(const char* seq, int seqlength, int k, int reach, int* results, int max_results) {
    int count = 0;
    for (int i = 0; i < seqlength - 2*k; i += k) {
        const char* seq1 = seq + i;
        if (contains_N_or_homogeneous(seq1, k)) continue;
        for (int j = i + k; j < i + reach && j < seqlength - k; j++) {
            const char* seq2 = seq + j;
            if (contains_N_or_homogeneous(seq2, k)) continue;
            int t = check_palindrome(seq1, seq2, k);
            if (t > 0) {
                if (count + 2 > max_results) return count; // 避免越界
                results[count*5 + 0] = i;
                results[count*5 + 1] = i + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = j;
                count++;
                results[count*5 + 0] = j;
                results[count*5 + 1] = j + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = i;
                count++;
            }
        }
    }
    return count;
}

// 分块扫描回文
int find_palindromes_chunk(const char* seq, int seqlength, int k, int reach, int start_idx, int end_idx, int* results) {
    int count = 0;
    if (start_idx < 0) start_idx = 0;
    if (end_idx > seqlength) end_idx = seqlength;
    for (int i = start_idx; i < end_idx - k; i += k) {
        const char* seq1 = seq + i;
        if (contains_N_or_homogeneous(seq1, k)) continue;
        for (int j = i + k; j < i + reach && j < end_idx - k; j++) {
            const char* seq2 = seq + j;
            if (contains_N_or_homogeneous(seq2, k)) continue;
            int t = check_palindrome(seq1, seq2, k);
            if (t > 0) {
                results[count*5 + 0] = i;
                results[count*5 + 1] = i + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = j;
                count++;
                results[count*5 + 0] = j;
                results[count*5 + 1] = j + k;
                results[count*5 + 2] = t;
                results[count*5 + 3] = k;
                results[count*5 + 4] = i;
                count++;
            }
        }
    }
    return count;
}