
// Last MOD: Carlo Sinigaglia, 2 Feb 2022
/* Tensors are defined A_{ijk} =  grad.phi_i grad.phi_j      phi_k
                       Bx_{ijk} = phi_i      gradx.phi_j     phi_k
                       By_{ijk} = phi_i      grady.phi_j     phi_k
                       C_{ijk} =  phi_i      phi_j           phi_k
*/
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

#ifdef _OPENMP
    #include <omp.h>
#else
    #warning "OpenMP not enabled. Compile with mex ADR_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);        /* Problem dimension 2,3 */
    int dim     = (int)(dim_ptr[0]);        
    int noe     = mxGetN(prhs[1]);             /* Number of elements, just get dimension */ 
    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[1]);
    int nln3    = nln*nln*nln;
    


    /* get operator to assemble */
    char *OP_string = mxArrayToString(prhs[8]);
    int OP[3] = {0, 0, 0};
    int    dT = 0;
    printf("Assembling %s operator \n", OP_string);
    if (strcmp(OP_string, "diffusion")==0)
    {
        OP[0] = 1;
    }
    
    if (strcmp(OP_string, "transport_x")==0)
    {
        OP[1]  = 1; /* dT=0 already defined above */ 
    }

    if (strcmp(OP_string, "transport_y")==0)
    {
        OP[1]  = 1;
        dT = 1; /* Transport in y direction */
    }
    
    if (strcmp(OP_string, "reaction")==0)
    {
        OP[2] = 1;
    }
    
    mxFree(OP_string);




    /*plhs pointers left-hand side*/
    plhs[0] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL);

    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAtens    = mxGetPr(plhs[2]);
    double* myAcoef    = mxGetPr(plhs[3]);
    
    
    int k,l,m;
 
    double C_t[dim];
    double C_d[dim][dim];

    for (k = 0; k < dim; k = k + 1 )
    {
        for (l = 0; l < dim; l = l + 1 )
        {
            C_d[k][l] = 0;
        }
        C_t[k] = 0;
    }
    

    for (l = 0; l < dim; l = l + 1 )
    {
        C_d[l][l] = 1;
    }

    
    for (l = 0; l < dim; l = l + 1 )
    {
        C_t[l] = 1;
    }

    int q;
    int NumQuadPoints     = mxGetN(prhs[3]);
    double* w   = mxGetPr(prhs[3]);
    double* invjac = mxGetPr(prhs[4]);
    double* detjac = mxGetPr(prhs[5]);
    double* phi = mxGetPr(prhs[6]);
    double* gradrefphi = mxGetPr(prhs[7]);

    /* does not exploit symmetry ? */

    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
                
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        

        int d1, d2;
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }



        int iii = 0;
        int ii = 0;
        int a, b,c;
    
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                for (c=0; c < nln; c = c+1)
                {

                    double aloc = 0;
                    for (q = 0; q < NumQuadPoints; q = q + 1 )
                    {
                        double diffusion = 0;
                        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                diffusion = diffusion + C_d[d1][d2] * phi[c+q*nln] * gradphi[d1][b][q] * gradphi[d2][a][q];
                            }
                        }

                        double transport = 0;
                        transport = transport + phi[c+q*nln] * gradphi[dT][b][q] * phi[a+q*nln];
                        
                        double reaction  = phi[c+q*nln] * phi[b+q*nln] * phi[a+q*nln];
                        
                        aloc = aloc + (OP[0] * diffusion + OP[1] * transport + OP[2] * reaction) * w[q];
                }

                    myArows[ie*nln3+iii] = elements[a+ie*numRowsElements];
                    myAcols[ie*nln3+iii] = elements[b+ie*numRowsElements];
                    myAtens[ie*nln3+iii] = elements[c+ie*numRowsElements];
                    myAcoef[ie*nln3+iii] = aloc*detjac[ie];
                    
                    iii = iii + 1;
                }
            }
                  
    
            ii = ii + 1;
        }
        
    }


}

