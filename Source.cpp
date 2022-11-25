#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<time.h>
#include <iostream>
int ProcNum;

int ProcRank;

int flag = 0;

int Size;

double *A; double *B; double* C;

using namespace std;

void PrintMatrix(double* pMatrix, int Size) {
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++)
			cout << pMatrix[i * Size + j] << " ";
		cout << endl;
	}
	cout << endl;
}
		


void RandInit(double* pMatrix, int Size) {
	//srand(100);
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) pMatrix[i * Size + j] = rand() / double(1000);
	}
}
void InitProcess(double*& A, double*& B, double*& C, int& S) {
	int Size=S;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0) {
		do {
			cout<<"\n--Square matrix multiplication--";
			cout<<"\nPlease, enter matrix size: "; Size = S++;
			if (Size < ProcNum) cout << "Matrix size is less than the number of processes! \n";
			if (Size % ProcNum != 0) cout << "Matrix size should be dividable by the number of processes! \n";
		} while ((Size < ProcNum) || (Size % ProcNum != 0));
	}
	S = Size;
	cout << "Size:" << Size << endl;
	if (Size < 10) flag = 1;
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		A = new double[Size * Size];
		B = new double[Size * Size];
		C = new double[Size * Size];
		RandInit(A, Size); RandInit(B, Size);
	}
}

void Transpose(double*& B, int dim) {
	double temp = 0;
	for (int i = 0; i < dim; i++)
		for (int j = i + 1; j < dim; j++) {
			temp = B[i * dim + j]; 
			B[i * dim + j] = B[j * dim + i]; 
			B[j * dim + i] = temp;
		}
}


void MatrixMultiplicationMPI(double*& A, double*& B, double*& C, int Size) {
	int dim = Size;
	int i, j, k, p, ind;
	double temp;
	MPI_Status Status;
	int ProcPartSize = dim / ProcNum;
	int ProcPartElem = ProcPartSize * dim;
	double* bufA = new double[dim * ProcPartSize];
	double* bufB = new double[dim * ProcPartSize];
	double* bufC = new double[dim * ProcPartSize];
	int ProcPart = dim / ProcNum, part = ProcPart * dim;
	if (ProcRank == 0) {
		Transpose(B, Size);
	}
	MPI_Scatter(A, part, MPI_DOUBLE, bufA, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, part, MPI_DOUBLE, bufB, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	srand(time(0));
	temp = 0.0;
	for (i = 0; i < ProcPartSize; i++)
	{
		for (j = 0; j < ProcPartSize; j++) {
			for (k = 0; k < dim; k++) temp += bufA[i * dim + k] * bufB[j * dim + k];
			bufC[i * dim + j + ProcPartSize * ProcRank] = temp; 
			temp = 0.0;
		}
	}
	int NextProc; int PrevProc;
	for (p = 1; p < ProcNum; p++) {
		NextProc = ProcRank + 1;
		if (ProcRank == ProcNum - 1) NextProc = 0;
		PrevProc = ProcRank - 1;
		if (ProcRank == 0) PrevProc = ProcNum - 1;
		MPI_Sendrecv_replace(bufB, part, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
		temp = 0.0;
		for (i = 0; i < ProcPartSize; i++) {
			for (j = 0; j < ProcPartSize; j++) {
				for (k = 0; k < dim; k++) {
					temp += bufA[i * dim + k] * bufB[j * dim + k];
				}
				if (ProcRank - p >= 0)
					ind = ProcRank - p;
				else ind = (ProcNum - p + ProcRank);
				bufC[i * dim + j + ind * ProcPartSize] = temp;
				temp = 0.0;
			}
		}
	}
	cout <<"Time: " << clock() / 1000.0 << endl;
	MPI_Gather(bufC, ProcPartElem, MPI_DOUBLE, C, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete[]bufA;
	delete[]bufB;
	delete[]bufC;
}

int main(int argc, char*argv[])
{
	MPI_Init(NULL, NULL);
	int s = atoi(argv[1]);
	InitProcess(A, B, C, s);
	MatrixMultiplicationMPI(A, B, C, s);
	MPI_Finalize();
	return 0;
}