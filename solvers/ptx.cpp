#ifdef USE_PASTIX6

#include <unistd.h>
#include <iostream>
#include <string>
#include <pastix.h>
#include <spm.h>
#include <cmath>

#include <mkl_spblas.h>
#include <mkl.h>
#include <omp.h>

extern "C" void ptx(void) {}

extern "C" void ptx_init(pastix_data_t **pastix_data_p,
                         spmatrix_t    **spm_p,
                         pastix_int_t   **iparm_p,
                         double         **dparm_p,
                         MPI_Fint       *comm_p)
{
    MPI_Comm pastix_comm = MPI_Comm_f2c(*comm_p);

    pastix_int_t    *iparm;
    double          *dparm;


    *iparm_p = new pastix_int_t[IPARM_SIZE];
    *dparm_p = new double[DPARM_SIZE];

    iparm = *iparm_p;
    dparm = *dparm_p;

    *spm_p = new spmatrix_t; // initialize pointer to data instance
    spmInit(*spm_p);

    pastixInitParam(iparm, dparm);
    
    dparm[DPARM_RELATIVE_ERROR] = 1e-8;
    dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
    iparm[IPARM_ITERMAX] = 40;
    dparm[DPARM_EPSILON_MAGN_CTRL] = 1e-24;
    
    int pastix_nthrd = 1;

#ifdef _OPENMP
#pragma omp parallel
  {
    #pragma omp master
    {
      pastix_nthrd = omp_get_num_threads();
    }
  }
#endif   

/*    
    if (pastix_nthrd>1){
      
      int procpernode, packsize, pe;
      MPI_Comm_rank(pastix_comm, &pe);
      
      std::cout<<"pastix_nthrd = "<<pastix_nthrd<<std::endl;
      iparm[IPARM_THREAD_NBR] = pastix_nthrd;
    
      int *bindtab = new int[iparm[IPARM_THREAD_NBR]];
      for(int i=0; i<iparm[IPARM_THREAD_NBR]; i++ ) {
        bindtab[i] = pe%pastix_nthrd + i;
      }
      pastixInitWithAffinity(pastix_data_p, pastix_comm, iparm, dparm, bindtab);
      delete [](bindtab);   
    }
    else{

      pastixInit(pastix_data_p, pastix_comm, iparm, dparm);
    }
*/
    if (pastix_nthrd>1){iparm[IPARM_THREAD_NBR] = pastix_nthrd;}
    pastixInit(pastix_data_p, pastix_comm, iparm, dparm);

    return;
}


extern "C" void ptx_set_mat(spmatrix_t **spm_, int *indx, int *n_p, int *nnz_p, int *nd_p, int *nnzd_p, int *dof_p,
                               int **rptr, int **cptr, double **values, int **loc2glob, int **glob2loc,
                               MPI_Fint *comm_, bool *update, bool *check)
{

    spmatrix_t *spm = *spm_;
    int n_d = *nd_p; int nnz_d = *nnzd_p;
    int gN = *n_p;  int gnnz = *nnz_p; int dof = *dof_p;
    int pe, npe;
    int rc = 0;
    
    MPI_Comm pastix_comm=MPI_Comm_f2c(*comm_);
    MPI_Comm_rank(pastix_comm, &pe);
    MPI_Comm_size(pastix_comm, &npe);    

    if (!(*update)){
      
      std::cout<<pe<<": PaStiX: setting new matrix, n_d, nnz_d "<<n_d<<" "<<nnz_d<<std::endl;

      spm->comm = pastix_comm;
      spm->baseval = *indx;
      spm->fmttype = SpmCSC;
      spm->clustnum = pe;
      spm->clustnbr = npe;

      spm->mtxtype = SpmGeneral;
      spm->flttype = SpmDouble;
      spm->n       = n_d;
      spm->nnz     = nnz_d;
      spm->layout  = SpmColMajor;
      spm->dof     = dof;

      spm->rowptr  = *rptr;
      spm->colptr  = *cptr;
      spm->values  = *values;

      spm->gN = gN;
      spm->gnnz = gnnz;
      spm->gNexp = gN;
      spm->gnnzexp = gnnz;
      spm->loc2glob = *loc2glob;
      spm->glob2loc = *glob2loc;

    }else{
      std::cout<<pe<<": PaStiX: updating matrix values "<<n_d<<" "<<nnz_d<<std::endl;
      spm->rowptr  = *rptr;
      spm->colptr  = *cptr;
      spm->values  = *values;
    }
    
    spmUpdateComputedFields(spm);
    
    //if (*check){
    //  spmatrix_t spm2;    
    //  rc = spmCheckAndCorrect(spm, &spm2);
    //  if (rc != 0) {        
    //    spmExit(spm);
    //    *spm = spm2;
    //    rc = 0;
    //  }
    //}

    return;
}

extern "C" void ptx_analyze(pastix_data_t **pastix_data_p, spmatrix_t **spm_p){

    spmatrix_t *spm = *spm_p;
    pastix_data_t *pastix_data = *pastix_data_p;

    int rc = 0;
    spmatrix_t spm2;
    rc = spmCheckAndCorrect(spm, &spm2);
    if (rc != 0) {
        spmExit(spm);
        spm = &spm2;
        rc = 0;
    }

    pastix_task_analyze(pastix_data, spm);
    return;
}


extern "C" void ptx_factorize(pastix_data_t **pastix_data_p, spmatrix_t **spm_p){

    spmatrix_t *spm = *spm_p;
    pastix_data_t *pastix_data = *pastix_data_p;

 // Perform the numerical factorization
    pastix_task_numfact(pastix_data, spm);
  
    return;

}

extern "C" void ptx_solve(pastix_data_t **pastix_data_p, spmatrix_t **spm_p, double **rhs_p, bool *refine){

    int nrhs = 1;
    void *b, *x, *x0 = NULL;
    int rc = 0;
    size_t size;

    double *rhs_d;

    spmatrix_t *spm = *spm_p;
    pastix_data_t *pastix_data = *pastix_data_p;

    // distribute RHS
    rhs_d = new double[spm->n];
    for (int i=0; i<spm->n; i++){
      rhs_d[i] = (*rhs_p)[spm->loc2glob[i] - spm->baseval];
    }

    size = pastix_size_of(spm->flttype)*spm->n*nrhs;
    x = malloc(size);
    memcpy(x, rhs_d, size);

    pastix_task_solve(pastix_data, nrhs, x, spm->n);
    
    if (*refine){
      pastix_task_refine(pastix_data, spm->n, nrhs, rhs_d, spm->n, x, spm->n);
    }

    //rc = spmCheckAxb(1e-6, nrhs, spm, x0, spm->n, rhs_d, spm->n, x, spm->n);    
   
    for (int i=0; i<(spm->gN); i++) {(*rhs_p)[i]=0.0;}; // prepare rhs to store the solution
#pragma omp for      
    for (int i=0; i<(spm->n); i++)
    {
      (*rhs_p)[spm->loc2glob[i] - spm->baseval] = ((double *)x)[i];
    }

    // the solution becomes available to all ranks
    MPI_Allreduce(MPI_IN_PLACE, (*rhs_p), spm->gN, MPI_DOUBLE, MPI_SUM, spm->comm);

    free(x);
    delete [](rhs_d);
    

    return;
}

extern "C" void ptx_finalize(pastix_data_t **pastix_data_p, spmatrix_t **spm_p,
                            pastix_int_t   **iparm_p,
                            double         **dparm_p){
  
    pastix_data_t *pastix_data = *pastix_data_p;
    spmatrix_t *spm = *spm_p;
  
    spm->colptr = NULL;
    spm->rowptr = NULL;
    spm->values = NULL;
    spm->loc2glob = NULL;
    spm->glob2loc = NULL;
    spm->dofs = NULL;  
   
    spmExit(spm);
    delete [](spm);
    delete [](*iparm_p);
    delete [](*dparm_p);
    //free(spm);
    //pastixFinalize(pastix_data_p);
    pastixFinalize(&pastix_data);
  
    return;
}

extern "C" void get_residual(int *n_p, int *nnz_p,
            int **irn_p, int **jcn_p,
            double **val_p, double **x_p,
            double **rhs_p, int *indx)
{
  double *val=*val_p, *x=*x_p, *rhs=*rhs_p;
  int *irn=*irn_p, *jcn=*jcn_p;
  int n=*n_p, nnz=*nnz_p, indexing = *indx;

  int  i;
  double rnorm=0.0, bnorm=0.0;
  double *res = new double[n];

  for (int i=0; i<n; i++){
    res[i] = 0.0;
  }

  // calculate matrix-vector product A*x
  for (int i=0; i<nnz; i++){
    res[irn[i]-indexing] += val[i]*x[jcn[i]-indexing];
  }

  for (int i=0; i<n; i++){
          bnorm+=rhs[i]*rhs[i];
          rnorm+=(rhs[i]-res[i])*(rhs[i]-res[i]);
  }

  delete [] res;

  std::cout<<"Relative error: "<<sqrt(rnorm/bnorm)<<std::endl;
  return;
}

#endif
