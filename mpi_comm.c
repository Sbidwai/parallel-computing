#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
typedef struct {float r; float i;} comp;
static comp temp;
#define swap(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define N 512
void fft(comp *r, int n, int is)
{
    int m,i,i1,j,k,i2,l,l1,l2;
    float c1,c2,z;
    comp t,u;
    if (is==0) return;
    i2 = n>>1;
    j=0;
    for(i=0;i<n-1;i++) {
        if(i < j)
        swap(r[i],r[j]);
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j+=k;
    }
    for (i=n,m=0;i>1; m++,i/=2);
    c1=-1.0;
    c2=0.0;
    l2=1;
    for (l=0;l<m;l++) {
        l1=l2;
        l2<<= 1;
        u.r= 1.0;
        u.i = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                i1 = i + l1;
                t.r = u.r * r[i1].r - u.i * r[i1].i;
                t.i = u.r * r[i1].i + u.i * r[i1].r;
                r[i1].r = r[i].r - t.r;
                r[i1].i = r[i].i - t.i;
                r[i].r += t.r;
                r[i].i += t.i;
            }
            z =  u.r * c1 - u.i * c2;
            u.i = u.r * c2 + u.i * c1;
            u.r = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (is == -1) 
        c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    if (is == 1) {     
        for (i=0;i<n;i++) {
            r[i].r/=n;
            r[i].i/=n;
        } }}
void getData(char file[15], comp **data){
    FILE *fp = fopen(file, "r");

    int i, j, res;

    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            res= fscanf(fp,"%g",&data[i][j].r);
            data[i][j].i=0.00;
        }  }
    fclose(fp);
}
void transp(comp **data, comp ** tr){
    int i, j;
    for (i=0;i<N; i++)
      for(j=0 ; j< N ; j++)
         tr[j][i] = data[i][j];
}
void m_point(comp **data1, comp **data2, comp **data3){
    int i, j;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            data3[i][j].r = (data1[i][j].r * data2[i][j].r) - (data1[i][j].i * data2[i][j].i);
            data3[i][j].i = (data1[i][j].r * data2[i][j].i) + (data1[i][j].i * data2[i][j].r);
        }
    }
}

void print(char file[15], comp **data){
    FILE *fp = fopen(file, "w");
    int i, j;
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++){
            fprintf(fp,"   %.7e",data[i][j].r);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}
int malloc2dfloat(comp ***array, int n, int m) {
    comp *p = (comp *)malloc(n*m*sizeof(comp));
    if (!p) return -1;
    (*array) = (comp **)malloc(n*sizeof(comp*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    int i;
    for (i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);
    return 0;
}

int free2dfloat(comp ***array) {
    free(&((*array)[0][0]));
    free(*array);
    return 0;
}

int main(int argc, char **argv){
    int my_rank, p;
    comp **data1, **data2, **data3, **data4;
    malloc2dfloat(&data1, N, N);
    malloc2dfloat(&data2, N, N);
    malloc2dfloat(&data3, N, N);
    malloc2dfloat(&data4, N, N);
    comp *v;

    char file1[15] = "im1";
    char file2[15] = "im2";
    char file3[15] = "mpi_comm_out";

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Datatype mystruct;
    int b_len[2] = { 1, 1 };
    MPI_Aint indices[2] = { 0, sizeof(float) };
    MPI_Datatype old_types[2] = { MPI_FLOAT, MPI_FLOAT };
    MPI_Type_struct( 2, b_len, indices, old_types, &mystruct );
    MPI_Type_commit( &mystruct );
    int i,j;
    double beg, end;
   int offset;
    int tag = 345;
    int rows = N/p;
    int lb = my_rank*rows;
    int hb = lb + rows;
    printf("%d have lb = %d and hb = %d\n", my_rank, lb, hb);
   if(my_rank == 0){
        getData(file1, data1);
        getData(file2, data2);
       printf("\nStarting clock.\n");
        beg = MPI_Wtime();
    }
    MPI_Bcast(&data1[0][0],N*N,mystruct,0,MPI_COMM_WORLD);
    MPI_Bcast(&data2[0][0],N*N,mystruct,0,MPI_COMM_WORLD);
  v = (comp *)malloc(N * sizeof(comp));
    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data1[i][j];
        }
        fft(v, N, -1);
        for (j=0;j<N;j++) {
            data1[i][j] = v[j];
        }
    }
    free(v);
    v = (comp *)malloc(N * sizeof(comp));
   for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data2[i][j];
        }
        fft(v, N, -1);
        for (j=0;j<N;j++) {
            data2[i][j] = v[j];
        }
    }
    free(v);
    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Recv(data1[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(data2[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
            }
        }
    }else{
        for(j = lb; j < hb; j++){
            MPI_Send(&data1[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
            MPI_Send(&data2[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
        }}
    if(my_rank == 0){
        transp(data1, data3);
        transp(data2, data4);
        for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Send(&data3[j][0], N, mystruct, i, tag, MPI_COMM_WORLD);
                MPI_Send(&data4[j][0], N, mystruct, i, tag, MPI_COMM_WORLD);
            }
        }
    }else{
        for(j = lb; j < hb; j++){
            MPI_Recv(data3[j], N, mystruct, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(data4[j], N, mystruct, 0, tag, MPI_COMM_WORLD, &status);
        }
    }
    v = (comp *)malloc(N * sizeof(comp));
    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data3[i][j];
        }
        fft(v, N, -1);
        for (j=0;j<N;j++) {
            data3[i][j] = v[j];
        }   }
    free(v);
    v = (comp *)malloc(N * sizeof(comp));
    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data4[i][j];
        }
        fft(v, N, -1);
        for (j=0;j<N;j++) {
            data4[i][j] = v[j];
        } }
    free(v);
    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Recv(data3[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(data4[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
            } }
    }else{
        for(j = lb; j < hb; j++){
            MPI_Send(&data3[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
            MPI_Send(&data4[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
        }}
    if(my_rank == 0){
        transp(data3,data1);
        transp(data4,data2);
        m_point(data1, data2, data3);
    }
    MPI_Bcast(&(data3[0][0]),N*N,mystruct,0,MPI_COMM_WORLD);
    v = (comp *)malloc(N * sizeof(comp));
    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data3[i][j];
        }
       fft(v, N, 1);
       for (j=0;j<N;j++) {
           data3[i][j] = v[j];
        } }
    free(v);
    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Recv(data3[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
            } }
    } else{
        for(j =lb;j<hb;j++){
            MPI_Send(&data3[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
        }}
    if(my_rank == 0){
        transp(data3,data4);
        for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Send(&data4[j][0], N, mystruct, i, tag, MPI_COMM_WORLD);
            }  }
    } else{
        for(j = lb; j < hb; j++){
            MPI_Recv(data4[j], N, mystruct, 0, tag, MPI_COMM_WORLD, &status);
        } }
    v = (comp *)malloc(N * sizeof(comp));
    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            v[j] = data4[i][j];
        }
        fft(v, N, 1);
        for (j=0;j<N;j++) {
            data4[i][j] = v[j];
        }}
    free(v);
    if(my_rank == 0){
       for(i=1;i<p;i++){
            offset=i*rows;
            for(j = offset; j < (offset+rows); j++){
                MPI_Recv(data4[j], N, mystruct, i, tag, MPI_COMM_WORLD, &status);
            } }
    } else{
        for(j = lb; j < hb; j++){
            MPI_Send(&data4[j][0], N, mystruct, 0, tag, MPI_COMM_WORLD);
        } }
    if(my_rank == 0){
        transp(data4,data3);
        end = MPI_Wtime();
        printf("\nElapsed time = %lf s.\n",(end - beg));
        printf("--------------------------------------------\n");
    }
    MPI_Finalize();
    if(my_rank == 0){
        print(file3, data3);
    }
    free2dfloat(&data1);
    free2dfloat(&data2);
    free2dfloat(&data3);
    free2dfloat(&data4);
    return 0; }
