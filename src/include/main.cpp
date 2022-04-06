#include <bits/stdc++.h>
using namespace std;

void printVector(vector<int> const &input)
{
    int n = input.size();
    for (int i = 0; i < n; i++) {
        cout << input[i]<<" ";
        if (i < n - 1) {
            cout << ", ";
        }
    }
    cout << " ";
}


int findMaxSumSubmatrix(vector<vector<int>> const &mat)
{
    // base case
    if (mat.size() == 0) {
        return 0;
    }
 
    // `M × N` matrix
    int M = mat.size();
    int N = mat[0].size();
 
    // `S[i][j]` stores the sum of submatrix formed by row 0 to `i-1`
    // and column 0 to `j-1`
    int S[M+1][N+1];
 
    // preprocess the matrix to fill `S`
    for (int i = 0; i <= M; i++)
    {
        for (int j = 0; j <= N; j++)
        {
            if (i == 0 || j == 0) {
                S[i][j] = 0;
            }
            else {
                S[i][j] = S[i-1][j] + S[i][j-1] - S[i-1][j-1] + mat[i-1][j-1];
            }
        }
    }
 
    int maxSum = INT_MIN;
    int rowStart, rowEnd, colStart, colEnd;
 
    // consider every submatrix formed by row `i` to `j`
    // and column `m` to `n`
    for (int i = 0; i < M; i++)
    {
        for (int j = i; j < M; j++)
        {
            for (int m = 0; m < N; m++)
            {
                for (int n = m; n < N; n++)
                {
                    // calculate the submatrix sum using `S[][]` in O(1) time
                    int submatrix_sum = S[j+1][n+1] - S[j+1][m] - S[i][n+1] + S[i][m];
 
                    // if the submatrix sum is more than the maximum found so far
                    if (submatrix_sum > maxSum)
                    {
                        maxSum = submatrix_sum;
                        rowStart = i;
                        rowEnd = j;
                        colStart = m;
                        colEnd = n;
                    }
                }
            }
        }
    }
 
    cout << "The maximum sum submatrix is\n\n";
    for (int i = rowStart; i <= rowEnd; i++) {
        vector<int> row;
        for (int j = colStart; j <= colEnd; j++) {
            row.push_back(mat[i][j]);
        }
        printVector(row);
    }
 
    return maxSum;
}


int maxSubArray(vector<int> A)
{
    int maxSum = A[0];
    int n = A.size();
    int dp[A.size()];
    dp[0] = A[0];
    for (int i = 1; i < n; i++)
    {
        maxSum = max(A[i], maxSum + A[i]); //we caompare the value of the current element and the value after adding the current value in the sum, this is because if we are considering the previous elements then we add those values otherwise the current value is taken.  
        dp[i] = max(dp[i - 1], maxSum); //we only update the maximum sum value going forward.
    }
    return dp[A.size() - 1];
}

int sum_answer(vector<int> &arr, int i, int &max, int sum, vector<int> &memo,int n)
{
    if (i >= n) // base case
    {
        return 0;
    }
    if (memo[i] != -1) // we are fetching the value from the memoized array
    {
        cout<<"checking memo array" << endl;
        return memo[i];
    }
    if (sum < 0) // if the sum gets less than 0 we update the value to 0 again
    {
        sum = 0; //this is the case when we are trying to get the best answer possible
    }

    if (sum > max)
    {
        max = sum; //updating the max value for the next recursion
    }

    sum += arr[i];
    sum_answer(arr, i + 1, max, sum, memo,n); // recursion for the elements of the array

    if (sum > max) // updating the maximum sum we encountered so far
    {
        max = sum;
    }
    memo[i] = max; // saving the value in the memoized array
    return max;    // we return the maximum value we encounter
}




int main(int argc, char *argv[])
{
    if (strcmp(argv[1],"1")==0) //TASK 1
    {
            int n; //driver function begins
            cout << "enter dimensions: ";
            cin >> n;
            int val;
            vector<int> arr;
            cout << "enter elements of the array: ";
            for (int j = 0; j < n; j++)
            {
                cin >> val;
                arr.push_back(val);
            } //driver function ends
            int max_sum = -2147483647, left = 0, right = 0;
            n=arr.size();
            for (int i = 0; i < n; i++) //this loops find the leftmost element of our subarray
            {
                for (int j = i; j < n; j++)
                {
                    int sum = 0;
                    for (int k = i; k <= j; k++)  //this loop finds the rightmost element of the subarray 
                    {
                        sum += arr[k];
                        if (sum > max_sum)
                        {
                            max_sum = sum;
                            left = i; //updating the leftmost element
                            right = j; //updating the right most element
                        }
                    }
                }
            }
        cout << left+1 << " " << right+1 <<" " <<max_sum<<endl; ////output
        return 0;
        }
    else if (strcmp(argv[1],"2")==0) //TASK 2
    {
        int n; //driver function begins
        cout << "enter dimensions: ";
        cin >> n;
        int val;
        vector<int> arr;
        cout << "enter elements of the array: ";
        for (int j = 0; j < n; j++)
        {
            cin >> val;
            arr.push_back(val);
        } //driver function ends
    int maxSum = arr[0], left = 0, right = 0;
    n=arr.size();
    for (int i = 0; i < n; i++)
    {
        int sum = 0;
        for (int j = i; j < n; j++)
        {
            sum += arr[j];
            if (maxSum < sum)
            {
                maxSum = sum;
                left = i; //updating the left most element of the subarray 
                right = j; //updating the right most element of the subarray
            }
        }
    }
    cout << left+1 << " " << right+1 <<" " <<maxSum<<endl; //output
    return 0;

    }
    else if (strcmp(argv[1],"3a")==0) //this is TASK 3A
    {
        int n;
        cout << "enter dimensions: ";
        cin >> n;
        int val;
        vector<int> nums;
        cout << "enter elements of the array: ";
        for (int j = 0; j < n; j++)
        {
            cin >> val;
            nums.push_back(val);
        }
        int max= 0;
        n = nums.size();

        vector<int> memo(n + 1, -1); // memoization array

        int answer = sum_answer(nums, 0,max,0,memo,n ); // here we expect the answer to be the maximum sum from the sub array.
        if (answer == 0) //this is for the condition all the elements of the given array are negative 
        {
           answer = *max_element(nums.begin(), nums.end()); //so we instead find the maximum element of the given array 
        }
    cout << answer << endl; //output 
    return 0;
    }
    else if (strcmp(argv[1],"3b")==0) //TASK 3B
    {
        int n; //driver function begins
        cout << "enter dimensions: ";
        cin >> n;
        int val;
        vector<int> nums;
         cout << "enter elements of the array: ";
        for (int j = 0; j < n; j++)
        {
            cin >> val;
            nums.push_back(val);
        } //driver function ends
        int T[nums.size()]; //array to save the value of the sum
        n = nums.size();
        T[0] = nums[0];
        int m = T[0];
        int left = 0;
        int left_temp=0;
        int right = 0;

        for (int i = 1; i < n; i++)
        {
            // if()

            T[i] = max(nums[i], nums[i] + T[i - 1]); //comparing either the maximum sum uptill now or the current value of the array element.
            
            if(T[i-1]>0){
                T[i] = nums[i] + T[i - 1];
            } else {
                T[i] = nums[i];
                left_temp = i;
            }

            if (T[i] > m){
                m = T[i];
                right = i;
                left = left_temp;
            }
                
        }
        cout<<left+1<< " "<< right+1<< " "<<m << endl; //output
        return 0;
    }
    else if (strcmp(argv[1],"4")==0)   //TASK 4
    {
        int n, m; //driver function begins
        cout << "enter dimensions: ";
        cin >> n >> m;
        int val;
        vector<vector<int>> mat;
        for (int i = 0; i < n; i++)
        {
            cout << "enter elements of the 2d array: ";
            vector<int> row;
            for (int j = 0; j < m; j++)
            {
                
                cin >> val;
                row.push_back(val);
            }
            mat.push_back(row);
        }//driver function ends
        int maxSum = -2147483647; //asssigning the maximum value to INT_MIN
    int leftr = 0, leftc = 0, rightr = 0, rightc = 0;

    for (int i = 0; i < n; ++i) //left most row value
    {
        for (int j = 0; j < m; j++) //left most column value
        { 
            for (int k = i; k < n; k++) //right most row value 
            {
                for (int l = j; l < m; l++) //right most column value
                {
                    int curSum = 0;
                    for (int row = i; row <= k; row++) //computing the sum in the next two for loops
                    {
                        for (int col = j; col <= l; col++)
                        {
                            curSum += mat[row][col];
                            if (curSum > maxSum)
                            {
                                //updating the values of maximum sum, leftmost row and column and rightmost row and column.
                                maxSum = curSum;
                                leftr = i;
                                leftc = j;
                                rightr = row;
                                rightc = col;
                            }
                        }
                    }
                }
            }
        }
    }
    cout << leftr+1 << " " << leftc+1 << " "<< rightr+1<< " "<<rightc+1<< " "<<maxSum<<endl; //output
    return 0;
    }
    else if (strcmp(argv[1],"5")==0)
    {
        int n, m; //driver function begins
        cout << "enter dimensions: ";
        cin >> n >> m;
        int val;
        vector<vector<int>> mat;
        for (int i = 0; i < n; i++)
        {
            cout << "enter elements of the 2d array: ";
            vector<int> row;
            for (int j = 0; j < m; j++)
            {
                cin >> val;
                row.push_back(val);
            }
            mat.push_back(row);
        }//driver function ends
        // int leftr = 0, leftc = 0, rightr = 0, rightc = 0;

        cout << "\nThe maximum sum is " << findMaxSumSubmatrix(mat);
        return 0; 
    }
    else if (strcmp(argv[1],"6")==0)
    {
        int n, m; //driver code begins
        cout << "enter dimensions: ";
        cin >> n >> m;
        int val;
        vector<vector<int>> mat;
        for (int i = 0; i < n; i++)
        {
            cout << "enter elements of the 2d array: ";
            vector<int> row;
            for (int j = 0; j < m; j++)
            {
                cin >> val;
                row.push_back(val);
            }
            mat.push_back(row);
        } //driver code ends
        int ar = INT_MIN; //assigning the value of the area to minimum integer 
        for (int l = 0; l < m; l++)
        {
            vector<int> sum(n);
            for (int r = l; r < m; r++)
            {
                for (int row = 0; row < n; row++)
                {
                    sum[row] += mat[r][row];
                }
                ar = max(ar, maxSubArray(sum));
            }
        }
        cout << ar<< endl; //output
        return 0;
    }
    else 
    cout << "The task does not exist" <<endl; 
    return 0;
}