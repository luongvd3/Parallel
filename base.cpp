// #include <mpi.h>
#include<iostream>
#include <stdio.h>
#include<bits/stdc++.h>
#include <sys/time.h>
using namespace std; 

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
} 

int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}
int max(int a, int b);
 
/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int lcs( std::string X, std::string Y, int m, int n )
{
int **L = new2d (m+1, n+1);
printf("%d %d\n",m,n);
int i, j;
 
/* Following steps build L[m+1][n+1] in bottom up fashion. Note
    that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
for (i=0; i<=m; i++)
{
    for (j=0; j<=n; j++)
    {
    if (i == 0 || j == 0)
        L[i][j] = 0;
 
    else if (X[i-1] == Y[j-1])
        L[i][j] = L[i-1][j-1] + 1;
 
    else
        L[i][j] = max(L[i-1][j], L[i][j-1]);
    }
}

/* L[m][n] contains length of LCS for X[0..n-1] and Y[0..m-1] */
return L[m][n];
}
 
/* Utility function to get max of 2 integers */
int max(int a, int b)
{
    return (a > b)? a : b;
}
 
/* Driver program to test above function */
int main()
{
// char X[] = "AGGTAB";
// char Y[] = "GXTXAYB";
std::string X, Y;
std::cin >> X;
std::cin >> Y;
int m = X.length(); 
int n = Y.length();
uint64_t start = GetTimeStamp(); 
printf("Length of LCS is %d\n", lcs( X, Y, m, n ) );
// print the time taken to do the computation
printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start)); 
return 0;
}