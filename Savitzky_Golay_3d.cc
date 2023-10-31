#include "Savitzky_Golay_3d.h"
#include <iostream>

using namespace std;

int Filter3d::alloc()
{
    filter =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    data = (double *) fftw_malloc(sizeof(double)*(size));
    datat = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    datat_filtered =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));

    plan_data2filter =   fftw_plan_dft_r2c_3d(z_size, y_size, x_size, data, filter, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_data2datat  =   fftw_plan_dft_r2c_3d(z_size, y_size, x_size, data,  datat, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_filtered2data = fftw_plan_dft_c2r_3d(z_size, y_size, x_size, datat_filtered, data,FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    return 0;
}

Filter3d::Filter3d(int nx, int ny, int nz,
                   int ld,//for derivatives
                   int m, //
                   int xsize, int ysize, int zsize):x_size(xsize),y_size(ysize), z_size(zsize)
{
    int m1 = m + 1;
    int x_len = nx * 2 + 1; //symmetric on both sides
    int y_len = ny * 2 + 1;
    int z_len = nz * 2 + 1;
    int col_len = x_len * y_len * z_len;
    size = x_size * y_size * z_size;
    size2 = (x_size/2 + 1) * y_size * z_size;
    alloc();
    box_range.init(x_len, y_len, z_len);
    filter_range.init(xsize, ysize, zsize);

    int m11 = m1*(m1+1)*(m1+2)/6;
    double *jacobian = new double[ m11 * col_len ];
    double *normal_matrix = new double[ m11 * m11];
    for(int u = -nz; u <=nz; u++){
       for(int t = -ny; t <= ny; t++){
           for(int s = -nx; s <= nx; s++){
               int ii = box_range(s, t, u);
               cout << "ii "<<ii<<endl;
               int col = 0;
               double ui = 1.;
               for (int i = 0; i < m1; i++){
                   double tj = 1.;
                   for(int j = 0; j < m1-i; j++){
                       double sk = 1.;
                       for(int k = 0; k < m1-i-j; k++, col++){
                           jacobian[ ii + col* col_len]
                           = sk* tj * ui;
                           sk*=s;
                       }
                       tj *=t;
                   }
                   ui *= u;
               }
            }
        }
    }
    for(int i = 0; i < m11 ; i++){
        for(int j = 0; j < col_len; j++){
            cout << jacobian[j+ i*col_len]<<" ";
        }cout << endl;
    }

	
	double double_one = 1.;
	double beta = 0.;
    int row_len = m11;
	cblas_dgemm(CblasColMajor,  CblasTrans,  CblasNoTrans, row_len, row_len, col_len, double_one, jacobian, col_len, jacobian, col_len,
				beta, normal_matrix, row_len); // Normal = J*J
	
	int *ipiv = new int[row_len]; 
	int info;
	dgetrf_(&row_len, &row_len, normal_matrix, &row_len, ipiv, &info); //rf normal
	double *work;
	double worksize;
	int lwork = -1;	
	dgetri_(&row_len, normal_matrix,&row_len, ipiv, &worksize, &lwork, &info);
	lwork = (int) (worksize);
	work = new double [lwork];
	dgetri_(&row_len, normal_matrix,&row_len, ipiv, work, &lwork, &info);
	double *cof=new double [row_len * col_len];
	cblas_dgemm(CblasColMajor,  CblasNoTrans, CblasTrans, row_len, col_len, row_len, double_one, normal_matrix, row_len, jacobian, col_len,
				beta, cof, row_len);
     

    for(int i = 0; i < size; i++) data[i] = 0;
    for(int k = -nz; k <= nz; k++){
        for(int j = -ny; j <= ny; j++){
            for(int i = -nx; i <= nx; i++) {
                cout <<box_range(i, j, k)<<" ";
                double temp = data[filter_range(i, j, k)] //box_range?
                = cof[row_len*(box_range(i, j, k))+ld];
                cout << temp <<" ";
            }
            cout << endl;
        }cout <<"------------------"<< endl;
    }
	fftw_execute(plan_data2filter);
}


int
Filter3d::apply(double *indata, double *outdata, int len)
{
  assert(len <= size);
  int i = 0;
  for(i = 0; i < len; i++){
     data[i] = indata[i];
  }
  while(i < size) { data[i++] = 0; }
  fftw_execute(plan_data2datat);

  for(int i =0; i < size2; i++) {
	datat_filtered[i][0] = (filter[i][0]*datat[i][0] - filter[i][1]*datat[i][1]) /size;
	datat_filtered[i][1] = (filter[i][0]*datat[i][1] + filter[i][1]*datat[i][0]) /size;
  }
  fftw_execute(plan_filtered2data); /* repeat as needed */

  for(i = 0; i < len; i++){
     outdata[i] = data[i];
  }
    return 0;
}
