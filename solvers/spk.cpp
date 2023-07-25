// Script for solving linear system Ax=B using STRUMPACK
// function spk can be called from the external code
// phase indicates the step (solver initialization,
// reordering and factorization, solving, and finalizing)
// I. Holod 29.01.2021
#ifdef USE_STRUMPACK
#include <iostream>
#include <string>
#include "hdf5.h"
#include <math.h>

#include "StrumpackSparseSolverMPIDist.hpp"
#include "sparse/CSRMatrix.hpp"
#include <chrono>

#include <sys/sysinfo.h>
#include <unistd.h>

#ifdef INTSIZE64
#define int_all int64_t
#else
#define int_all int
#endif

using namespace strumpack;

extern "C" void spk(void) {}

extern "C" void spk_init(StrumpackSparseSolverMPIDist<double,int_all>** spss_, int** iparm_, MPI_Fint* comm_) {

  StrumpackSparseSolverMPIDist<double,int_all>* spss= *spss_;
  int* iparm = *iparm_;
  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int thread_level,rank,P;
  double eps=1e-36, epsr=1.e-12;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);
  MPI_Query_thread(&thread_level);
  if (thread_level < MPI_THREAD_FUNNELED && rank == 0)
    std::cout << "MPI implementation does not support MPI_THREAD_FUNNELED"
              << std::endl;

  *spss_= new StrumpackSparseSolverMPIDist<double,int_all>(comm);
  spss = *spss_;

  if (iparm[0] == 1) {
    spss->options().set_Krylov_solver(KrylovSolver::DIRECT);
  } else if (iparm[0] == 2) {
    spss->options().set_Krylov_solver(KrylovSolver::PREC_GMRES);
  } else if (iparm[0] == 3) {
    spss->options().set_Krylov_solver(KrylovSolver::REFINE);
  }
  
  if (iparm[1] == 1) {
    spss->options().set_reordering_method(ReorderingStrategy::METIS);
  } else if (iparm[1] == 2) {
    spss->options().set_reordering_method(ReorderingStrategy::PARMETIS);
  } else if (iparm[1] == 3) {
    spss->options().set_reordering_method(ReorderingStrategy::PTSCOTCH);
  }
  // set reordering to METIS and MatchingJob to MAX_DIAGONAL_PRODUCT_SCALING as needed for older version
  if (STRUMPACK_VERSION_MAJOR<7) {iparm[1] = 1; iparm[2] = 1;}
  spss->options().enable_METIS_NodeND();

  if (iparm[2] == 0) {
    spss->options().set_matching(MatchingJob::NONE);
  } else if (iparm[2] == 1) {
    spss->options().set_matching(MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING);
  }

  spss->options().set_rel_tol(epsr);
  spss->options().set_abs_tol(eps);
  spss->options().set_maxit(200);
  spss->options().set_gmres_restart(50);
  spss->options().set_verbose(true);

//  spss->options().set_compression(CompressionType::HSS);
//  spss->options().set_compression_min_sep_size(512);
//  spss->options().set_compression_min_front_size(1024);
//  spss->options().HSS_options().set_rel_tol(1e-6);
//  spss->options().HSS_options().set_abs_tol(1e-10);
//  spss->options().BLR_options().set_rel_tol(1e-4);
//  spss->options().BLR_options().set_abs_tol(1e-8);

  return;
}


extern "C" void spk_set_mat(int_all* n_, int_all** dist_, int_all** irn_, int_all** jcn_, double** val_,
        StrumpackSparseSolverMPIDist<double,int_all>** spss_,MPI_Fint* comm_,bool* upd_) {
// set and factorize (distributed) matrix

  int_all n=*n_;
  int_all *dist=*dist_;
  int_all *irn=*irn_;
  int_all *jcn=*jcn_;
  double *val=*val_;
  bool upd=*upd_;

  StrumpackSparseSolverMPIDist<double,int_all>* spss= *spss_;

  std::chrono::steady_clock::time_point t0, t1;
  t0 = std::chrono::steady_clock::now();
#ifdef NEWSPK
  if (upd){
    std::cout<<"Updating matrix values"<<std::endl;
    bool symmetric_pattern = true;
    spss->update_matrix_values(n, irn, jcn, val, dist, symmetric_pattern);
  }else
#endif
  {
    spss->set_distributed_csr_matrix(n, irn, jcn, val, dist);
  }

  t1 = std::chrono::steady_clock::now();
  std::cout<<"Time to set matrix (s) = "<< std::chrono::duration_cast<
       std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}


extern "C" void spk_reord(StrumpackSparseSolverMPIDist<double,int_all>** spss_,MPI_Fint* comm_) {
// reorder (distributed) matrix

  StrumpackSparseSolverMPIDist<double,int_all>* spss= *spss_;

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int rank,P;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);
  std::chrono::steady_clock::time_point t0, t1;

  // Reordering
  t0 = std::chrono::steady_clock::now();
  spss->reorder();
  t1 = std::chrono::steady_clock::now();
  if (!rank)
    std::cout<<"Time to reorder (s) = "<< std::chrono::duration_cast<
    std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}


extern "C" void spk_fact(StrumpackSparseSolverMPIDist<double,int_all>** spss_,MPI_Fint* comm_) {
// factorize (distributed) matrix

  StrumpackSparseSolverMPIDist<double,int_all>* spss= *spss_;

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int rank,P;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);
  std::chrono::steady_clock::time_point t0, t1;

  // Factorization
  t0 = std::chrono::steady_clock::now();
  spss->factor();
  t1 = std::chrono::steady_clock::now();
  if (!rank)
      std::cout<<"Time to factorize (s) = "<< std::chrono::duration_cast<
      std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;

  return;
}


extern "C" void spk_solve(int_all* n_, int_all ** dist_, double** rhs_,
        StrumpackSparseSolverMPIDist<double,int_all>** spss_,MPI_Fint* comm_,int* phase) {

  int_all n=*n_;
  double* rhs=*rhs_;
  int_all *dist = *dist_;

  StrumpackSparseSolverMPIDist<double,int_all>* spss= *spss_;

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int thread_level,rank,P;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);
  std::chrono::steady_clock::time_point t0, t1;

  int_all n_local = dist[rank+1]-dist[rank];

  // set local RHS
  std::vector<double> b(n_local), x(n_local);

#pragma omp for
  for (int_all i=dist[rank]; i<dist[rank+1]; i++)
    b[i-dist[rank]]=rhs[i];

  t0 = std::chrono::steady_clock::now();
  spss->solve(b.data(),x.data(),false);

// Gather the solution
  std::vector<double> x_glob(n), x_buf(n);
  x_glob.assign(n,0);
  x_buf.assign(n,0);

#pragma omp for
  for (int_all i=dist[rank]; i<dist[rank+1]; i++)
    x_buf[i]=x[i-dist[rank]];

  MPI_Allreduce(x_buf.data(), x_glob.data(), n, MPI_DOUBLE_PRECISION, MPI_SUM, comm);

  t1 = std::chrono::steady_clock::now();
  if (!rank){
    std::cout<<"Time to solve (s) = "<< std::chrono::duration_cast<
      std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;
    }

#pragma omp for
  for (int_all i=0;i<n;i++){
    (*rhs_)[i] = x_glob[i];
  }

  x.clear();
  b.clear();
  x_glob.clear();
  x_buf.clear();

  return;
}

extern "C" void spk_delete_factors(StrumpackSparseSolverMPIDist<double,int>** spss_) {
  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;

  spss->delete_factors();

  return;
}

extern "C" void spk_finalize(StrumpackSparseSolverMPIDist<double,int_all>** spss_,MPI_Fint* comm_) {

  delete *spss_;
  scalapack::Cblacs_exit(1);
  return;
}


int_all* distribute(int_all n, int P){
  int_all* dist;
  int* nl;

  dist = new int_all[P+1];
  nl = new int[P];

  for (int i=0; i<P; i++){
    nl[i] = floor(n/P);
    if (i<n%P)
      nl[i]+=1;
  }
  dist[0]=0;
  for (int i=0; i<P; i++){
    dist[i+1] = dist[i] + nl[i];
  }
  return dist;
}


#endif
