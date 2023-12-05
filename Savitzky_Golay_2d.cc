#include "Savitzky_Golay_2d.h"
#include <iostream>

using namespace std;

int Filter2d::alloc()
{
    filter =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    data = (double *) fftw_malloc(sizeof(double)*(size));
    datat = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    datat_filtered =(fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * size2));
    plan_data2filter = fftw_plan_dft_r2c_2d(y_size, x_size, data, filter, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_data2datat =  fftw_plan_dft_r2c_2d(y_size, x_size, data, datat, FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    plan_filtered2data = fftw_plan_dft_c2r_2d(y_size, x_size, datat_filtered, data,FFTW_ESTIMATE);//the last dimension has the fastest-varying index in the array
    return 0;
}

void
Filter2d::init_polynom(int in_nx, int in_ny, int in_ld, int in_m){
    nx = in_nx;
    ny = in_ny;
    ld = in_ld;
    m = in_m;
    init_polynom();
 //   init_filter();
}

void
Filter2d::init_filter()
{
    for(int i = 0; i < size; i++) data[i] = 0;
    for(int j = -ny; j <= ny; j++){
        for(int i = -nx; i <= nx; i++) {
            double temp = data[filter_range(i, j)] //box_range?
            = cof[m11*(box_range(i, j))+ld];
        }
    }
    fftw_execute(plan_data2filter);
}

void
Filter2d::init_filter(int xsize, int ysize)
{
    x_size = xsize;
    y_size = ysize;
    size = x_size * y_size ;
    size2 = (x_size/2 + 1) * y_size;
    alloc();
    filter_range.init(xsize, ysize);
    init_filter();
}

void 
Filter2d::init_polynom()
{
    int m1 = m+1;
    int x_len = nx*2 + 1;
    int y_len = ny*2 + 1;
    col_len = x_len * y_len;
    box_range.init(x_len, y_len);
    m11 = m1*(m1+1)/2;
    double *jacobian = new double[ m11 * col_len ];
    double *normal_matrix = new double[ m11 * m11];
    for(int k = -ny; k <= ny; k++){
        for(int s = -nx; s <= nx; s++){
            int ii = box_range(s, k);
            int col = 0;
            for (int i = 0; i < m1; i++){
                for(int j = 0; j < m1-i; j++, col++){
                    double temp = jacobian[ ii + col* col_len]
                    = pow (k, i) * pow(s, j);
                }
            }
        }
    }
    cout << "jacobian"<<endl;
    int mat_size = col_len*col_len;
    double double_one = 1.;
    double beta = 0.;
    int row_len = m11;
    cblas_dgemm(CblasColMajor,  CblasTrans,  CblasNoTrans, row_len, row_len, col_len, double_one, jacobian, col_len, jacobian, col_len,
                beta, normal_matrix, row_len); // Normal = J*J
    
    cout << "normal"<<endl;
    int *ipiv = new int[row_len];
    int info;
    dgetrf_(&row_len, &row_len, normal_matrix, &row_len, ipiv, &info); //rf normal
    double worksize;
    int lwork = -1;
    dgetri_(&row_len, normal_matrix,&row_len, ipiv, &worksize, &lwork, &info);
    lwork = (int) (worksize);
    double *work = new double [lwork];
    dgetri_(&row_len, normal_matrix,&row_len, ipiv, work, &lwork, &info);
    if(cof) delete[] cof;
    cof=new double [row_len * col_len];
    cblas_dgemm(CblasColMajor,  CblasNoTrans, CblasTrans, row_len, col_len, row_len, double_one, normal_matrix, row_len, jacobian, col_len,
                beta, cof, row_len);
    cout <<"deleting"<<endl;
    delete [] work;
    delete [] ipiv;
    delete [] normal_matrix;
    delete [] jacobian;
}


Filter2d::Filter2d(int nx, int ny, int ld, int m, int xsize, int ysize):x_size(xsize),y_size(ysize)
{
//    double *temp = reinterpret_cast<double *>( fftw_malloc(sizeof(double) * size));

    cout << "constructing"<<endl;
    init_polynom(nx, ny, ld, m);
    cout << " polynom OK"<<endl;
    init_filter(xsize, ysize);
    cout << "init filter "<<endl;
}


int
Filter2d::apply(double *indata, double *outdata, int len)
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
