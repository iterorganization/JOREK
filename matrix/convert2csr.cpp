#include <iostream>
#include <string>
#include <math.h>
#include <chrono>
#include <vector>
#include<algorithm>
#ifdef USEMKL
#include "mkl_spblas.h"
#endif

#ifdef INTSIZE64
#define int_all int64_t
#else
#define int_all int
#endif

#if (defined(USEMKL)) && (!defined(INTSIZE64))
extern "C" void convert2csr(int *indx, int *n_, int *m_, int *nnz_, int **irn, int **jcn, double **val)
{
  int nr = *n_, mc = *m_, nnz = *nnz_;
  std::cout<<"nr = "<<nr<<" mc = "<<mc<<" nnz = "<<nnz<<std::endl;

  sparse_index_base_t indexing = sparse_index_base_t(*indx);
  sparse_matrix_t cooA;
  sparse_matrix_t csrA;  
  sparse_status_t stat;


  std::chrono::steady_clock::time_point t0, t1;
  t0 = std::chrono::steady_clock::now();

// create mkl coordinate sparse matrix
  mkl_sparse_d_create_coo(&cooA, indexing, nr, mc, nnz, *irn, *jcn, *val);  

// convert to csr format  
  mkl_sparse_convert_csr(cooA, SPARSE_OPERATION_NON_TRANSPOSE, &csrA);
  mkl_sparse_destroy(cooA);
  mkl_sparse_order(csrA); // important

// export csr values, rowptr (Begin and End counting) and colind
  MKL_INT *rowptrE, *rowptrB, *colind;
  double *values;
  mkl_sparse_d_export_csr(csrA, &indexing, &nr, &mc, &rowptrB, &rowptrE, &colind, &values);
  
  nnz = rowptrE[nr-1] - (*indx); 

  if (nnz!=(*nnz_)) 
	  std::cout<<"New nnz: "<<nnz<<" Old nnz "<< *nnz_<<std::endl;    
  *nnz_ = nnz;

#pragma omp for  
  for (int i=0; i<nr; i++){
	  (*irn)[i] = rowptrB[i];
  }
  
  (*irn)[nr] = nnz + (*indx);

#pragma omp for
  for (int i=0; i<nnz; i++){
	  (*jcn)[i] = colind[i];
	  (*val)[i] = values[i];
  }

  mkl_sparse_destroy(csrA);

  t1 = std::chrono::steady_clock::now();
  std::cout<<"coo2csr (MKL) (s) = "<< std::chrono::duration_cast<
             std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}

#else

extern "C" void convert2csr(int *indx_, int_all *n_, int_all *m_, int_all *nnz_, int_all **irn, int_all **jcn, double **val)
{
  int indx=*indx_;
  int_all nc =*m_, nr=*n_, nnz=*nnz_;

  std::vector<int_all> rptr(nr+1 ,0);

  std::cout<<"rows = "<<nr<<" cols = "<<nc<<" nnz = "<<nnz<<std::endl;
  if (nr>nnz+1){std::cout<<"Error: #rows > nnz+1"<<std::endl; return;}

  std::chrono::steady_clock::time_point t0, t1;
  t0 = std::chrono::steady_clock::now();

  for (int_all i=0; i<nnz; i++){
    rptr[(*irn)[i] + 1 - indx] += 1;
  }

  (*irn)[0] = indx;
  (*irn)[nr] = nnz + indx;
  for (int_all i=1; i<nr; i++){(*irn)[i] = rptr[i] + (*irn)[i-1];}
  rptr.clear();

  //for (int_all i=0; i<nnz; i++){(*jcn)[i] -= indx;}
    
  //*indx_ = 0;

  t1 = std::chrono::steady_clock::now();
  std::cout<<"coo2csr (no-MKL) (s) = "<< std::chrono::duration_cast<
			std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}
#endif

extern "C" void sortunique(int_all *nnz_, int64_t **ijn)
{
  int_all nnz=*nnz_;
  std::chrono::steady_clock::time_point t0, t1;

  t0 = std::chrono::steady_clock::now();
  
  std::vector<int64_t> vidx((*ijn),(*ijn)+nnz);
  
  std::sort(vidx.begin(), vidx.end());
  
  vidx.erase(std::unique(vidx.begin(), vidx.end()), vidx.end());
  nnz = vidx.size();

  *nnz_ = nnz;
  for (int_all i=0; i<nnz;i++){(*ijn)[i] = vidx[i];}
  vidx.clear();

  t1 = std::chrono::steady_clock::now();
  std::cout<<"Sorting time (s) = "<< std::chrono::duration_cast<
             std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}

