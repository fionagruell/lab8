// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void b(double& b1,double& b2,double& b3,double& b4,double Theta);
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,double* k1,double* k2,double* k3,double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.3,x=0;
	const double L = 20*M_PI;
    double y0[dim],yi[dim];
    
	double yn[dim];
	double b1,b2,b3,b4;
	double k1[dim], k2[dim], k3[dim], k4[dim];
	
 	for(double j=0;j<5;j+=0.001){
	 y0[0]=(1./10) +j;
	 y0[1]=0;
	x=0;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx,k1,k2,k3,k4);
          	if(yn[1]<0 && y0[1]>0) break;
                for(int i=0; i<dim; i++) y0[i] = yn[i];
		//out << x << "\t" << y0[0] << "\t" << y0[1] <<  endl;
	}
	double ThetaL=0,ThetaR=1,Theta=0.5;
	yi[1]=1;
       while(abs(yi[1])>1e-6)
	{
	  b(b1,b2,b3,b4,Theta);
	  //yi[0]=y0[0]+dx*(b1*k1[0]+b2*k2[0]+b3*k3[0]+b4*k4[0]);
	  yi[1]=y0[1]+dx*(b1*k1[1]+b2*k2[1]+b3*k3[1]+b4*k4[1]);
	  if(yi[1]>0)
	    ThetaL=Theta;
	  else ThetaR=Theta;
	  Theta=(ThetaL+ThetaR)/2.;
	  //cout << Theta << "\t" << yi[1] << endl;
	} 
	
	//cout << Theta << endl;
 out <<(1./10) +j <<"\t"<<x+Theta*dx -dx <<"\t"<<y0[0]<<"\t"<<y0[1]<<endl;
	}
	out.close();
	return(0);
	}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,double* k1,double* k2,double* k3,double* k4)
{
	const int dim = 2;
	

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	double y[2] = { y0[0], y0[1] };

  y0[0] = y[1];
  y0[1] = -y[0]/(sqrt(1+pow(y[0],2)));
	
}
void b(double& b1,double& b2,double& b3,double& b4,double Theta){
b1= Theta-3*Theta*Theta/2.+2.*Theta*Theta*Theta/3.;
b2=Theta*Theta-2.*Theta*Theta*Theta/3.;
b3=b2;
b4=-Theta*Theta/2.+2.*Theta*Theta*Theta/3.;
  
}