/********************************************
The code solves a 2D Barkley model equation in a square domain using finite difference method
with the parameters:

          epsilon = 0.02
          a       = 0.75
          b       = 0.01
          L       = 20,   domain
          tfinal = 40
          
          outputs: u and v

Written by: Yao Gahounzo
      
      03/13/2021
******************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

// *** Function interfaces ***

// getArgs_mpi reads input parameters Nfrom command line
void getArgs(int *N, int argc, char *argv[]);

// write results to file
void write_results(double *x, double *y, double *u, int N);

//declaration od functions f and g
double funf(double ep, double a, double b,double u, double v);
double gunf(double u, double v);

// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, idn;
  int id_lft, id_rgt, id_top, id_bot;

  double eps = 0.02; 
  double a = 0.75;
  double b = 0.01;
  double f,g;
   
  //read command line arguments
  getArgs(&N, argc, argv);

  // define domain size (assume both x and y directions are the same)
  double L = 150.0;
  double tfinal = 40.0;
  
  //compute the grid spacing
  double h = 2*L/(N-1);
  // compute the time step
  double dt = 0.0001;
  int ntime = (int)tfinal/dt + 1;
  
  //allocate arrays
  double *u       = malloc((N+2)*(N+2)*sizeof(double));
  double *v       = malloc((N+2)*(N+2)*sizeof(double));
  double *u_new   = malloc((N+2)*(N+2)*sizeof(double));
  double *v_new   = malloc((N+2)*(N+2)*sizeof(double));

  double *x       = malloc((N+2)*sizeof(double));
  double *y       = malloc((N+2)*sizeof(double));

  //initialize x and y
  for(i=1;i<N+1;i++){
    x[i] = - L + (i-1)*h;
    y[i] = - L + (i-1)*h;
  }
  
  //initialize array u
  for(j=1;j<N+1;j++){
     for(i=1;i<N+1;i++){
    	id = ID_2D(i,j,N);

    	if(y[j] >= 0){

    	    u[id] = 1;
    	}else if(y[j] < 0){
    	    u[id] = 0;
    	}

    	if(x[i] >= 0){

    	    v[id] = 1;
    	}else if(x[i] < 0){
    	    v[id] = 0;
    	}
			
      }
  }
  
  // time integration

  for(int n = 1; n < ntime; n++){

    //initialize BC in ghost cells

    for(i=1;i<N+1;i++){
    id = ID_2D(i,0,N); //bottom ghost
    idn = ID_2D(i,1,N);
    u[id] = u[idn];
    
    id = ID_2D(i,N+1,N); //top ghost
    idn = ID_2D(i,N,N);
    u[id] = u[idn];

    id = ID_2D(0,i,N); //left ghost
    idn = ID_2D(1,i,N);
    u[id] = u[idn];

    id = ID_2D(N+1,i,N); //right ghost
    idn = ID_2D(N,i,N);
    u[id] = u[idn];
    }
    
    //go over interior points 2:N-1
    for(j=1;j < N+1;j++){
      for(i=1;i < N+1;i++){
  	    id = ID_2D(i,j,N);
  	    id_lft = ID_2D(i-1,j,N); //index left
  	    id_rgt = ID_2D(i+1,j,N); //index righ
  	    id_bot = ID_2D(i,j-1,N); //index bottom
  	    id_top = ID_2D(i,j+1,N); //index top
	    
  	    f = funf( eps, a, b, u[id], v[id]);
        u_new[id] = u[id] + dt*f + (dt/(h*h))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

        g = gunf(u[id], v[id]);
        v_new[id] = v[id] +  dt*g;
      }
    }
    
    //enforce Dirichlet BC by copying values from ghosts
    for(i=1;i<N+1;i++){
      j = 1;
      id = ID_2D(i,j,N);
      id_lft = ID_2D(i-1,j,N); //index left
      id_rgt = ID_2D(i+1,j,N); //index righ
      id_bot = ID_2D(i,j-1,N); //index bottom
      id_top = ID_2D(i,j+1,N); //index top
      
      f = funf( eps, a, b, u[id], v[id]);
      u_new[id] = u[id] + dt*f + (dt/(h*h))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

      g = gunf(u[id], v[id]);
      v_new[id] = v[id] +  dt*g;
      
    } 

    for(i = 1; i < N+1; i++){
      j = N+1;
      id = ID_2D(i,j,N);
      id_lft = ID_2D(i-1,j,N); //index left
      id_rgt = ID_2D(i+1,j,N); //index righ
      id_bot = ID_2D(i,j-1,N); //index bottom
      id_top = ID_2D(i,j+1,N); //index top
      
      f = funf( eps, a, b, u[id], v[id]);
      u_new[id] = u[id] + dt*f + (dt/(h*h))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

      g = gunf(u[id], v[id]);
      v_new[id] = v[id] +  dt*g;
    }

    for(j = 1; j < N+1; j++){
      i = 1;
      id = ID_2D(i,j,N);
      id_lft = ID_2D(i-1,j,N); //index left
      id_rgt = ID_2D(i+1,j,N); //index righ
      id_bot = ID_2D(i,j-1,N); //index bottom
      id_top = ID_2D(i,j+1,N); //index top
      
      f = funf( eps, a, b, u[id], v[id]);
      u_new[id] = u[id] + dt*f + (dt/(h*h))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

      g = gunf(u[id], v[id]);
      v_new[id] = v[id] +  dt*g;
    }

    for(j = 1; j < N+1; j++){
      i = N;
      id = ID_2D(i,j,N);
      id_lft = ID_2D(i-1,j,N); //index left
      id_rgt = ID_2D(i+1,j,N); //index righ
      id_bot = ID_2D(i,j-1,N); //index bottom
      id_top = ID_2D(i,j+1,N); //index top

      f = funf( eps, a, b, u[id], v[id]);
      u_new[id] = u[id] + dt*f + (dt/(h*h))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

      g = gunf(u[id], v[id]);
      v_new[id] = v[id] +  dt*g;
    }
    
    //update solution 
    for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
	      id = ID_2D(i,j,N);
	      u[id] = u_new[id];
	      v[id] = v_new[id];
	
      }
    }
      
  } // end time integration

  write_results(x, y, u, N);
  
    
#ifdef DEBUG // print the content of array u including ghosts
  for(j=0; j<N+2; j++){
    for(i=0; i<N+2; i++){
      id = ID_2D(i,j,N);
      printf("%f\t",u[id]);
    }
    printf("\n");
  }
#endif
  
  free(u);
  free(v);
  free(u_new);
  free(v_new);
  free(x);
  free(y);
 
  return 0;
}


// *** Function definitions ***
void getArgs(int *N, int argc, char *argv[])
{

    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./Barkley N \n");
       
       }
    else
      {
	//Use input argument
	*N = atoi(argv[1]);	
      }
}

// function f(u,v)

double funf(double ep, double a, double b,double u, double v){
  

  double f = (1.0/ep)*u*(1.0 - u)*(u - ((v + b)/a));

  return f;
}

// function g(u,v)

double gunf(double u, double v){
  

  double g = u - v;

  return g;
}

// write to a file
void write_results(double *x, double *y, double *u, int N){
  int i,j, id;

  FILE *f;
  f = fopen("barkley_serial.dat","wt"); //open file
  
  if ( !f ){
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "WRITE_RESULTS - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }

  for(j=1; j<N+1; j++){
    for(i=1; i<N+1; i++){
      id = ID_2D(i,j,N);
      fprintf(f,"   %24.16g\t%24.16g\t%24.16g\n",x[i],y[j],u[id]);
    }
  }
  fclose(f);
  
}
