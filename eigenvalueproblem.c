#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "lapacke.h"
#include <math.h>
#include <time.h>



int maximum = 351;
float h = 0.1;

//Potential energy = 1/2k(x-x_0)
void count(float *array, float *v,int maximum) //Create function and locate to a already defined array 
{
   
   FILE *f=fopen("pos_0.1_test.txt", "w");
   FILE *fd=fopen("pot_0.1_test.txt", "w");
    for (int i = 0;i < maximum; i++)
    {
        array[i] = (i*h)-round(maximum/2)*h; 
        v[i] = array[i]*array[i]+200*exp(-(array[i]*array[i])*2);
        fprintf(f,"%f\n",array[i]);
        fprintf(fd,"%f\n",v[i]);

    }
        fclose(f); 
        fclose(fd); 
}


void hamil(float hamiltonian[maximum][maximum], int maximum, float potent[maximum], float A[maximum*maximum], float hw, float c_00, float c_n11, float c_n22, float c_11, float c_22) //Create function and locate to a already defined array 
{
    int shift = 0;
    for (int k = 0; k < maximum; k++){
        for (int j = 0; j < maximum; j++){
            hamiltonian[k][j] = 0;
            A[k*(maximum-1)+j] = 0;
    
        }
        
    }
    for (int k = 0; k<maximum*maximum; k ++ )
    {
        A[k]  = 0;
    }
    hamiltonian[0][0] = ((c_00*1) + potent[0])* hw+ shift;
    hamiltonian[0][1] = (c_11*1 )* hw;
    hamiltonian[0][2] = (c_22*1)* hw;
    A[0] = (c_00*1 + potent[0])* hw + shift;
    A[1] = (c_11* hw);
    A[2] = (c_22* hw);

    hamiltonian[1][0] = c_n11 * hw;
    hamiltonian[1][1] = ((c_00) + potent[1])* hw + shift;
    hamiltonian[1][2] = (c_11*1)* hw;
    hamiltonian[1][3] = (c_22*1 )* hw;
    A[maximum] = (c_n11)* hw;
    A[maximum +1] = ((c_00*1) + potent[1])* hw+ shift;
    A[maximum +2] = (c_11)* hw;
    A[maximum +3] = (c_22)* hw;

    hamiltonian[maximum-2][maximum-2] = ((c_00*1)+ potent[maximum-2])* hw+ shift;
    hamiltonian[maximum-2][maximum-1] =  (c_11*1)* hw;
    hamiltonian[maximum-2][maximum-3] =  (c_n11*1)* hw;
    hamiltonian[maximum-2][maximum-4] =  (c_n22*1)* hw;
    A[maximum*(maximum-1)-2] = ((c_00*1)+ potent[maximum-2])* hw+ shift;
    A[maximum*(maximum-1)-1] = (c_11*1)* hw;
    A[maximum*(maximum-1)-3] = (c_n11*1)* hw;
    A[maximum*(maximum-1)-4] = (c_n22*1)* hw;

    hamiltonian[maximum-1][maximum-1] = ((c_00*1+ potent[maximum-1]))* hw+ shift;
    hamiltonian[maximum-1][maximum-2] =  (c_n11*1)* hw;
    hamiltonian[maximum-1][maximum-3] =  (c_n22*1)* hw;
    A[maximum*maximum-1] = ((c_00*1)+ potent[maximum-1])* hw+ shift;
    A[maximum*maximum-2] = (c_n11*1)* hw;
    A[maximum*maximum-3] = (c_n22*1)* hw;
    for (int i = 2;i < maximum-2; i++){
        hamiltonian[i][i-2] = (c_n22*1)* hw; 
        hamiltonian[i][i-1] = (c_n11*1)* hw; 
        hamiltonian[i][i+1] = (1*c_11)* hw; 
        hamiltonian[i][i+2] = (1*c_22)* hw; 
        hamiltonian[i][i] = ((c_00*1) + potent[i])* hw+ shift;
        A[(i)*maximum + (i-2)]=  (c_n22*1)* hw;
        A[(i)*maximum + (i-1)]=  (c_n11*1)* hw;
        A[(i)*maximum + (i)] =  ((c_00*1) + potent[i])* hw+ shift;
        A[(i)*maximum + (i+1)] =(c_11)* hw; 
        A[(i)*maximum + (i+2)] =(c_22)* hw; 
        
    }
    A[maximum*maximum-maximum] = 0;
}


int main() 
{   

    float potential[maximum];
    int element[maximum];
    float arr[maximum]; // Define array
    count(arr, potential,maximum); //call function which will output the array
    float c_n2 = -1/(-12*h*h);
    float c_n1= 16/(-12*h*h);
    float c_0 = -30/(-12*h*h);
    float c_1 = 16/(-12*h*h);
    float c_2 = -1/(-12*h*h);    
    float hamilton[maximum][maximum];
    float A[maximum*maximum]; 
    float hw_2= 0.50;
    hamil(hamilton, maximum, potential,A,hw_2, c_0, c_n1, c_n2, c_1, c_2);
    //ssyev_2stage() tridiagonal
    //ssygv_2stage() solves VPsi = EPsi
/*
    float w[maximum];
    int lwork= 100*maximum; 
    float work[lwork];
    int info;
    int dim = maximum;
    int lda = maximum;
    ///DIAGONALIZE WITH LAPACK
    float A_unchanged[maximum*maximum];

    for (int i = 0; i <maximum*maximum;i++)
    {
        A_unchanged[i] = A[i];
    }
    clock_t begin = clock();

    ssyev_("V","U", &dim, A, &lda ,w, work, &lwork,&info);
    clock_t end= clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f\n", time_spent);
    FILE *g=fopen("eigenvalue_0.1_test.txt", "w");
    for (int i = 0;i<maximum; i++)
    {
        printf("%d",info);
        fprintf(g,"%f\n",w[i]);
        printf("\n");
        

    }  
    fclose(g); 
    
//Eigenvectors 
    FILE *e=fopen("eigenfunctions_0.1_test.txt", "w");
    for (int i = 0;i<maximum*maximum; i++)
    {
        
        printf("\n");
        fprintf(e,"%f\n",A[i]);
        
        }
    fclose(e);
    
    //check that we get correct eigenvectors 
    float left_side[maximum];
    for (int i = 0; i <maximum;i++)

    {
        for (int j=0; j<maximum;j++)
        {   
            left_side[i]=A_unchanged[(i*maximum)+j]*A[j]+left_side[i];
        }
        printf(" left side: %f and right side: %f\n",left_side[(i)],A[i]*w[0]);
    }
*/

    //////////////Power iteration//////////////
    
    float y_i[maximum];
    float y_im1[maximum];

    for (int i = 0; i<maximum;i++)
    {
       y_i[i] = rand()%10+1;

    } 
    int order = maximum;
    int nrhs = 1;
    int ipiv;
    int lwork = maximum;
    //int lwork = maximum * maximum + 5 * maximum;
    float work[lwork];
    int lda = maximum;
    int ldb = maximum;
    float A_unchanged[maximum * maximum];
    int rows = maximum;
    int coloumns = maximum;
    int JPVT[maximum];
    float RCOND = 0.00001;
    int RANK;
    int LWORK = maximum  * maximum;
    float WORK[LWORK];
    int INFO;
    float absolute = 0;
    float num = 0;  
    float den = 0;   
    float eigenvalue;

    
    clock_t begin = clock();
    for (int iteraive = 0; iteraive<10; iteraive=iteraive+1)
    {
        absolute = 0;
        RANK = maximum;
        for (int k = 0; k < maximum; k++) {     // norm the current eigenstate y_i
            absolute += y_i[k] * y_i[k];
            JPVT[k] = 0;
        }
        absolute = sqrt(absolute);
        for (int k = 0; k < maximum; k++) {
            y_i[k] = y_i[k] / absolute;
        }
        for (int i = 0; i < maximum * maximum; i++)
        {
            A_unchanged[i] = A[i];
        }       
        sgelsy_(&rows, &coloumns, &nrhs, A_unchanged, &lda, y_i, &ldb, JPVT, &RCOND, &RANK, WORK, &LWORK, &INFO);
       // printf("info : %d \n", INFO);
    }

    absolute = 0;
    RANK = maximum;
    for (int k = 0; k < maximum; k++) {     // norm the current eigenstate y_i
        absolute += y_i[k] * y_i[k];
        JPVT[k] = 0;
    }
    absolute = sqrt(absolute);
    for (int k = 0; k < maximum; k++) {
        y_i[k] = y_i[k] / absolute;
        y_im1[k] = y_i[k];          // copy y_i in y_{i-1}, iterate y_i once more and than calc eigenvalue
    }
    for (int i = 0; i < maximum * maximum; i++)
    {
        A_unchanged[i] = A[i];
    }
    sgelsy_(&rows, &coloumns, &nrhs, A_unchanged, &lda, y_i, &ldb, JPVT, &RCOND, &RANK, WORK, &LWORK, &INFO);
    for (int k = 0; k < maximum; k++) {
        num = y_i[k] * y_im1[k]+ num;
        den = y_im1[k] * y_im1[k]+ den;
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent:%f\n", time_spent);
    eigenvalue = den / num;
    printf("Eigenvalue: %f\n", eigenvalue);

    FILE *ef=fopen("ef_iterative.txt", "w");
    for (int i = 0;i<maximum*maximum; i++)
    {
        printf("\n");
        fprintf(ef,"%f\n",A[i]);
        
        }
    fclose(ef);



}
