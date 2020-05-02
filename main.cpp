// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Junda Lou

#include <iostream>
#include <vector>
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <string>
#include "mpi.h"

using namespace std;

void mmser(vector<vector<int> >, vector<vector<int> >);
void mm1d(vector<vector<int> >, vector<vector<int> >);
void mm2d(vector<vector<int> >, vector<vector<int> >);
vector<vector<int> > init_matrix(int, int, bool);
int get_random_number(int);
string matrix_to_string(vector<vector<int> >);

int main(int argc, char** argv) {
    int m, n, q;
    if (argc == 4) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        q = atoi(argv[3]);
    }

    srand(time(NULL));
    MPI_Init(NULL, NULL);
    
    int threads;
    MPI_Comm_size(MPI_COMM_WORLD, &threads);
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // printf("P: %d\nMatrix Dimensions: (%d, %d), (%d, %d)\n", threads, m, n, n, q);

    vector<vector<int> > m1 = init_matrix(m, n, true);
    vector<vector<int> > m2 = init_matrix(n, q, true);
    // vector<vector<int> > m1{ {1, 2, 3}, {4, 5, 6} };
    // vector<vector<int> > m2{ {7, 8}, {9, 10}, {11, 12} };
    // vector<vector<int> > m1{ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
    // vector<vector<int> > m2{ {10, 11, 12}, {13, 14, 15}, {16, 17, 18} };

    mmser(m1, m2);
    mm1d(m1, m2);
    mm2d(m1, m2);

    // cout << "left\n" << matrix_to_string(m1) << "right\n" << matrix_to_string(m2) << "result\n" << matrix_to_string(mult);

    MPI_Finalize();
}

void mmser(vector<vector<int> > m1, vector<vector<int> > m2) {
    int m = m1.size(), n = m1[0].size(), q = m2[0].size();

    int thread_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);
    if (thread_id == 0) {
        vector<vector<int>> res = init_matrix(m, q, false);

        double timeSpent = 0.0;
        clock_t begin = clock();
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < q; j++) {
                    res[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        clock_t end = clock();
        timeSpent += (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Serial Elapsed Time: %.8f seconds\n", timeSpent);
    }
}

void mm1d(vector<vector<int> > m1, vector<vector<int> > m2) {
    int m = m1.size(), n = m1[0].size(), q = m2[0].size();

    int thread_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);

    if (thread_id < max(m, q)) {
        if (thread_id < m) {
            for (int j = 0; j < q; j++) {
                if (j != thread_id) {
                    MPI_Send(m1[thread_id].data(), n, MPI_INT, j, 0, MPI_COMM_WORLD);
                }
            }
        }
        
        if (thread_id < q) {
            for (int i = 0; i < m; i++) {
                int buff[n]; 
                if (i != thread_id) {
                    MPI_Recv(&buff, n, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int res_buff[3];
                    res_buff[0] = 0;
                    res_buff[1] = i;
                    res_buff[2] = thread_id;
                    for (int j = 0; j < n; j++) {
                        // res[i][thread_id] = buff[j] * m2[j][thread_id];
                        res_buff[0] += buff[j] * m2[j][thread_id];
                    }
                    MPI_Send(&res_buff, 3, MPI_INT, max(m, q), 0, MPI_COMM_WORLD);
                } else {
                    int res_buff[3];
                    res_buff[0] = 0;
                    res_buff[1] = thread_id;
                    res_buff[2] = thread_id;
                    for (int j = 0; j < n; j++) {  
                        // res[thread_id][thread_id] = m1[thread_id][j] * m2[j][thread_id];
                        res_buff[0] += m1[thread_id][j] * m2[j][thread_id];    
                    }
                    MPI_Send(&res_buff, 3, MPI_INT, max(m, q), 0, MPI_COMM_WORLD);
                }
            }
        }
    } else if (thread_id == max(m, q)) {
        double timeSpent = 0.0;
        clock_t begin = clock();

        vector<vector<int>> res = init_matrix(m, q, false);
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < m; j++) {
                int buff[3];
                MPI_Recv(buff, 3, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                res[buff[1]][buff[2]] = buff[0];
            }
        }

        clock_t end = clock();
        timeSpent += (double)(end - begin) / CLOCKS_PER_SEC;
        printf("MM1D Elapsed Time: %.8f seconds\n", timeSpent);
        // cout << matrix_to_string(res);
    }
}

void mm2d(vector<vector<int> > m1, vector<vector<int> > m2) {
    int m = m1.size(), n = m1[0].size(), q = m2[0].size();

    int thread_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);
    int tr = thread_id / n, tc = thread_id % n;
    if (thread_id < max(m, q) * n) {
        if (tr < m && tc < n) {
            for (int i = 0; i < q; i++) {
                if (thread_id != tc * q + i) {
                    MPI_Send(&m1[tr][tc], 1, MPI_INT, tc * q + i, 0, MPI_COMM_WORLD);
                }
            }
        }

        tr = thread_id / q, tc = thread_id % q;
        if (tr < n && tc < q) {
            for (int i = 0; i < m; i++) {
                if (thread_id != i * n + tr) {
                    int buff;
                    MPI_Recv(&buff, 1, MPI_INT, i * n + tr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int res_buff[3];
                    res_buff[0] = buff * m2[tr][tc];
                    res_buff[1] = i;
                    res_buff[2] = tc;
                    MPI_Send(&res_buff, 3, MPI_INT, max(m, q) * n, 0, MPI_COMM_WORLD);
                } else {
                    int res_buff[3];
                    res_buff[0] = m1[tr][tc] * m2[tr][tc];
                    res_buff[1] = tr;
                    res_buff[2] = tc;
                    MPI_Send(&res_buff, 3, MPI_INT, max(m, q) * n, 0, MPI_COMM_WORLD);
                }
            }
        }
    } else if (thread_id == max(m, q) * n) {
        double timeSpent = 0.0;
        clock_t begin = clock();

        vector<vector<int>> res = init_matrix(m, q, false);
        for (int i = 0; i < n * q; i++) {
            for (int j = 0; j < m; j++) {
                int buff[3];
                MPI_Recv(&buff, 3, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                res[buff[1]][buff[2]] += buff[0];
            }
        }
        clock_t end = clock();
        timeSpent += (double)(end - begin) / CLOCKS_PER_SEC;
        printf("MM2D Elapsed Time: %.8f seconds\n", timeSpent);
        // cout << matrix_to_string(res);
    }
}

vector<vector<int> > init_matrix(int row, int col, bool random) {
    vector<vector<int> > matrix;
    for (int i = 0; i < row; i++) {
        vector<int> tmp;
        for (int j = 0; j < col; j++) {
            tmp.push_back(random ? get_random_number(10) : 0);
        }
        matrix.push_back(tmp);
    }
    return matrix;
}

int get_random_number(int bound) {
    return rand() % bound + 1;
}

string matrix_to_string(vector<vector<int> > matrix) {
    string res = "";
    int m = matrix.size(), n = matrix[0].size();
    for (int i = 0; i < m; i++) {
        res += "[";
        for (int j = 0; j < n - 1; j++) {
            res += to_string(matrix[i][j]) + ", ";
        }
        if (n > 0) {
            res += to_string(matrix[i][n - 1]);
        }
        res += "]\n";
    }
    res += "\n";
    return res;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
