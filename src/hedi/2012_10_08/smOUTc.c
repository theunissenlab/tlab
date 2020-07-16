#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "sm.h"





int ODEfunc (double t, const double y[], double f[],void *params)
{

  ODEparams * pa = (ODEparams *)params;
  double alpha=gsl_spline_eval(pa->a_spline,t, pa->a_acc);
  double beta=gsl_spline_eval(pa->b_spline,t, pa->b_acc);
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];

  f[0] = yy;
  f[1] = alpha*gamma2 - beta*gamma2*x-gamma2*x*x*x - gamma*x*x*yy+gamma2*x*x - gamma*x*yy;
  return GSL_SUCCESS;
}

int ODEjac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{

  ODEparams * pa = (ODEparams *)params;
  double alpha=gsl_spline_eval(pa->a_spline,t, pa->a_acc);
  double beta=gsl_spline_eval(pa->b_spline,t, pa->b_acc);
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];


  //double * p = (double *)params;
  	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  	gsl_matrix * m = &dfdy_mat.matrix; 
	/*	0              1
		-k -2*c*y    -(p-b)-cx^2 
	*/
  	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
  	gsl_matrix_set (m, 1, 0, -beta*gamma2 - 3*gamma2*x*x - 2*gamma*x*yy + 2*gamma2*x - gamma*yy);
  	gsl_matrix_set (m, 1, 1, -gamma*x*x - gamma*x);
	dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  return GSL_SUCCESS;
}




void ODErun(double * data, double tmax, double dt,ODEparams * pa)
{

	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
  	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-10, 0.0);
  	gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2);
  	gsl_odeiv_system sys = {ODEfunc, ODEjac, 2, pa};
  	int totalframe=ceil(tmax/dt);
  	int i=0;
	double t1=dt,t=0.0;  		
	double h = 1e-10;
	double y[2] = { pa->x0,pa->y0};
	

	while(t1<tmax)
    {
       while (t < t1)
  	    {  
	               int status = gsl_odeiv_evolve_apply (e, c, s,
					       &sys, 
					       &t, t1,
					       &h, y);
	  
        	  if (status != GSL_SUCCESS)
	               break;
      	}	  
	    i++;	
      	t1+=dt;
      	data[i]=y[0];
//      	printf("%lf %lf\n",t1,y[0]);

	}
	gsl_odeiv_evolve_free (e);
  	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

}


void ODEfilter(double * data, int nsample,int nfft,ODEparams * pa)
{
	int i,k;
	gsl_fft_real_wavetable * real;
  	gsl_fft_halfcomplex_wavetable * hc;
  	gsl_fft_real_workspace * work;
  	work = gsl_fft_real_workspace_alloc (nfft);
	real = gsl_fft_real_wavetable_alloc (nfft);
  	hc = gsl_fft_halfcomplex_wavetable_alloc (nfft);

  	double * filter = (double *) calloc(nfft,sizeof(double));

	double * vdata = (double *) calloc(nsample+nfft/2,sizeof(double));
	memcpy(vdata,data,nsample*sizeof(double));

  	for (i = 0; i < nfft; i++)
    {
    	double ff=i*44100.0/nfft/2; 
		double fil=gsl_spline_eval(pa->f_spline,ff, pa->f_acc);
  		filter[i]=fil;
    }

    for(i=0;i<nsample;i+=nfft)
    {
    	for(k=0;k<nfft;k++)
    		{
    			(data+i)[k] *= (1-k/(1.*nfft))*k/(1.*nfft); /* hamming window*/
    			(vdata+i+nfft/2)[k] *= (1-k/(1.*nfft))*k/(1.*nfft);
			}
	  	gsl_fft_real_transform (data+i, 1, nfft,real, work);
    	gsl_fft_real_transform (vdata+i+nfft/2, 1, nfft,real, work);
    	for(k=0;k<nfft;k++)
    	{
    		(data+i)[k]=(data+i)[k]*filter[k];
    		(vdata+i+nfft/2)[k]=(vdata+i+nfft/2)[k]*filter[k];
    	}
    	gsl_fft_halfcomplex_inverse (data+i, 1, nfft, hc, work);
    	gsl_fft_halfcomplex_inverse (vdata+i+nfft/2, 1, nfft, hc, work);
    	for (k=0;k<nfft;k++)
    	 (data+i)[k]=((data+i)[k]+(vdata+i)[k])/2;
    }

}

void save_wave_file(char * fname,double * data,int size)
{
  SNDFILE * file; 
  SF_INFO sfinfo ;

  sfinfo.samplerate       = 44100 ;
  sfinfo.frames           = size ;
  sfinfo.channels         = 1 ;
  sfinfo.format           = (SF_FORMAT_WAV | SF_FORMAT_PCM_16) ;

  file=sf_open(fname,SFM_WRITE,&sfinfo);
  sf_write_double (file, data, sfinfo.channels * size);
  sf_close(file);
}

double * load_wave_file(char * fname, int * ns)
{
	SNDFILE * file; 
  	SF_INFO sfinfo ;
	file=sf_open(fname,SFM_READ,&sfinfo);
	(*ns)=sfinfo.frames;
	double * data=(double*)malloc((*ns)*sizeof(double));
	sf_readf_double(file,data,(*ns));
	sf_close(file);
	return data;
}

double * get_spectrum_data(double * data, int ns, int nfft, int ntarget, int nslices,int * chunks)
{

 // compute the spectrum that will be the target
	int nchunks=ns/nfft;
	if (ns%nfft!=0) nchunks += -1;
	*chunks=nchunks;
	double * spectrum = (double *)calloc(nslices*nchunks*ntarget,sizeof(double));
	memcpy(spectrum,data,ns*sizeof(double));
	double buffer[nfft];

	int i,offset=0;
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;
  	
  	work = gsl_fft_real_workspace_alloc (nfft);
	real = gsl_fft_real_wavetable_alloc (nfft);
  	for(i=0;i<nchunks*nfft;i+=nfft/nslices)
    {
    	memcpy(buffer,data+i,nfft*sizeof(double));
      	gsl_fft_real_transform (buffer, 1, nfft,real, work);
      	memcpy(spectrum+offset,buffer,ntarget*sizeof(double));
      	offset+=ntarget;
    }
    return spectrum;
}

ODEparams * ODEinit(int nap,double * alpha,double * talpha, int nbp, double* beta, double *tbeta,int nfp,double * filter, double * freqs,double gamma, double x0,double y0)
{
	ODEparams * pa= (ODEparams*)malloc(sizeof(ODEparams));

	pa->gamma=gamma;
	pa->x0=x0;
	pa->y0=y0;
	
	pa->nap=nap;
	pa->nbp=nbp;
	pa->nfp=nfp;

	pa->alpha=calloc(nap,sizeof(double));
	pa->talpha=calloc(nap,sizeof(double));
	pa->beta=calloc(nbp,sizeof(double));
	pa->tbeta=calloc(nbp,sizeof(double));
	pa->filter=calloc(nfp,sizeof(double));
	pa->freqs=calloc(nfp,sizeof(double));
	
	int i;
	memcpy(pa->alpha,alpha,nap*sizeof(double));
	memcpy(pa->talpha,talpha,nap*sizeof(double));
	memcpy(pa->beta,beta,nbp*sizeof(double));
	memcpy(pa->tbeta,tbeta,nbp*sizeof(double));
	memcpy(pa->filter,filter,nfp*sizeof(double));
	memcpy(pa->freqs,freqs,nfp*sizeof(double));
//	int i;
//	for(i=0;i<nbp;i++) printf("%lg",pa->tb[i]);
	pa->a_acc=gsl_interp_accel_alloc ();
	pa->b_acc=gsl_interp_accel_alloc ();
	pa->f_acc=gsl_interp_accel_alloc ();
	pa->a_spline = gsl_spline_alloc (gsl_interp_cspline,nap);
	gsl_spline_init (pa->a_spline, pa->talpha, pa->alpha, pa->nap);
	pa->b_spline = gsl_spline_alloc (gsl_interp_cspline,nbp);
	gsl_spline_init (pa->b_spline, pa->tbeta, pa->beta, pa->nbp);
	pa->f_spline = gsl_spline_alloc (gsl_interp_cspline,nfp);
	gsl_spline_init (pa->f_spline, pa->freqs, pa->filter, pa->nfp);

	pa->gamma=gamma;
	pa->x0=x0;
	pa->y0=y0;
	return pa;
}

void ODEupdate(ODEparams* pa,int nap,double * alpha, int nbp, double* beta,int nfp,double * filter)
{

	memcpy(pa->alpha,alpha,nap*sizeof(double));
	memcpy(pa->beta,beta,nbp*sizeof(double));	
	memcpy(pa->filter,filter,nfp*sizeof(double));
	gsl_spline_init (pa->a_spline, pa->talpha, pa->alpha, pa->nap);	
	gsl_spline_init (pa->b_spline, pa->tbeta, pa->beta, pa->nbp);
	gsl_spline_init (pa->f_spline, pa->freqs, pa->filter, pa->nfp);

}

void OPTresponse(double *p, double *x, int m, int n, void *data)
{

	OPTparams  * OPTpa=(OPTparams*)data;
	ODEparams * ODEpa=OPTpa->ODEpa;

	double * spectrum= OPTpa->spectrum;
	double * talpha = OPTpa->talpha;
	double * tbeta= OPTpa->tbeta;
	double * freqs= OPTpa->freq;
	
	int nap=11;
	int nbp=4;	
	int nfp=10;
	double * alpha=p;
	double * beta=p+nap;
	double * filter=p+nap+nbp;	
	ODEparams * pa=ODEupdate(ODEpa,nap,alpha,nbp,beta,nfp,filter);
	ODErun(x,tmax,dt,pa);
	ODEfilter(x,nsample,nfft,pa);

}

OPTparams * OPTinit(char * fname,int nap,double * salpha,int nbp,double sbeta, int nfp,double *sfreq)
{
	OPTparams * pa=(OPTparams*)malloc(sizeof(OPTparams));
	sprintf(pa->fname,"%s",fname);	
	pa->nfft=512;
	pa->samplerate=44100.;
	pa->cut_freq=8000.;

	pa->ntarget=floor(cut_freq/samplerate*nfft);
	pa->nslices=4;

	double * target_data=load_wave_file(fname, &pa->target_ns);
	double tmax=pa->target_ns/pa->samplerate;

	spectrum=get_spectrum_data(target_data,pa->target_ns,pa->nfft, pa->ntarget,pa->nslices,&pa->nchunks);

	pa->n=pa->nchunksn*pa->nslices*pa->ntarget;

	pa->nap=nap;
	pa->talpha=(double*)calloc(nap,sizeof(double));
	for(i=0;i<nap;i++)
		pa->talpha[i]=i*tmax/(nap-1);

	pa->nbp=nbp;
	pa->tbeta=(double*)calloc(nbp,sizeof(double));
	for(i=0;i<nbp;i++)
		pa->tbeta[i]=i*tmax/(nbp-1);

	pa->nfp=nfp;
	pa->freqs=(double*)calloc(nfp,sizeof(double));


	pa->ODEpa=ODEinit(nap,salpha,talpha,nbp,sbeta,tbeta,nfp,sfilter,freqs,double gamma, double x0,double y0))
}



int main(int argc, char ** argv)
{

	int target_ns;
	char fname[100];
	int i,j;
	sprintf(fname,"%s",argv[1]);
	double * target_data=load_wave_file(fname, &target_ns);


	int nfft=512,nchunks;
	double samplerate=44100.;
	double cut_freq=8000.;
	int ntarget=floor(cut_freq/samplerate*nfft);
	int nslices=4;

	double * spectrum=get_spectrum_data(target_data,target_ns,nfft, ntarget,nslices,&nchunks);
	for(i=0;i<nchunks*nslices;i++)
	{
		for(j=0;j<ntarget;j++)
			printf("%e ",spectrum[i*ntarget+j]);
		printf("\n");
	}
	int nap=11;
	int nbp=4;	
	int nfp=10;
	

	double tmax= 0.200;
	double gamma=24000;

	double alpha[]={-0.1129,-0.1474,-0.1542,-0.1552,-0.1522,-0.1507,-0.1486,-0.1447,-0.1293,-0.1291,-0.1187};
	double talpha[]={0.0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20};
	double beta[]={0.0105,0.0105,0.0105,0.0105};
	double tbeta[]={0.0,0.10, 0.15,0.20};
	double filter[]={10,100,1000,2000,8000,4000,1000,1.0, 0.0, 0.0};
	double freqs[]={0.0, 1000.0, 2000.0,3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0,samplerate};  	
	double dt=1/(samplerate);

	int totalframe=ceil(tmax/dt);

	ODEparams * pa=init(nap,alpha,talpha,nbp,beta,tbeta,nfp,filter,freqs,gamma,0,1);

	int nsample=totalframe;

	if (nsample%nfft!=0) nsample=nsample + nfft- (nsample%nfft);
	

  	double * data = (double*) calloc(totalframe,sizeof(double));
	runOde(data,tmax,dt,pa);
	runFilter(data,nsample,nfft,pa);
	
	// for(i=0;i<nsample;i++)
	// 	printf("%e\n",data[i]);


	return 0;
}