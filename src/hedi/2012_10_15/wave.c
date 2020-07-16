#include "wave.h"





void
save_wave_file(char *fname, double *data, int size)
{
	SNDFILE        *file;
	SF_INFO         sfinfo;

	sfinfo.samplerate = 44100;
	sfinfo.frames = size;
	sfinfo.channels = 1;
	sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_PCM_16);
	double          Mx = getABSmax(data, size);
	scaleData(1 / Mx, data, size);
	file = sf_open(fname, SFM_WRITE, &sfinfo);
	sf_write_double(file, data, sfinfo.channels * size);
	sf_close(file);
}

double         *
load_wave_file(char *fname, int *ns)
{
	SNDFILE        *file;
	SF_INFO         sfinfo;
	file = sf_open(fname, SFM_READ, &sfinfo);
	(*ns) = sfinfo.frames;
	double         *data = (double *) malloc((*ns) * sizeof(double));
	sf_readf_double(file, data, (*ns));
	sf_close(file);
	return data;
}

// typedef struct 
// {
// 	double * data; // whole data (including surroundings)

// } wave_slice;

// typedef struct 
// {
// 	wave_slice ** slices;
//  double * wave_data;
// 	int nslices; 
// 	int nfft;
// } wave_array;



wave_slice * create_slice(double * wave_data,int start,int nfft)
{
	wave_slice * ws = (wave_slice*)malloc(sizeof(wave_slice));
	ws->data = (double*)calloc(2*nfft,sizeof(double));
	memcpy(ws->data,wave_data+start,2*nfft*sizeof(double));
	int j=0; double m=0;
	for(j=0;j<2*nfft;j++) m+= ws->data[j]/(2*nfft);
	for(j=0;j<2*nfft;j++) ws->data[j] -= m;
	return ws;

}
wave_array * create_array(char * fname,int nfft)
{
	wave_array * wa=(wave_array*)malloc(sizeof(wave_array));
	double * data=load_wave_file(fname,&wa->ns);
	wa->nfft=nfft;
	// we assume 44100.0 1 channel
	// detecting the numbers of slices
	int nsli=wa->ns/nfft;
	if (wa->ns%nfft!=0) nsli++;
	wa->nslices=nsli;

	int ndata=nsli*nfft;
	wa->wave_data=(double *)calloc(ndata+nfft,sizeof(double));
	memcpy(wa->wave_data+nfft/2,data,wa->ns*sizeof(double));
	free(data);
	 double m =getABSmax(wa->wave_data, ndata+nfft);
//	 printf("MM %e\n",m);
	 scaleData(1/m,wa->wave_data,ndata+nfft);
	// taking care of each slices 
	int i;
	wa->slices=(wave_slice**)malloc(nsli*sizeof(wave_slice*));
	for(i=0;i<nsli;i++)
	{
		int start=i*nfft;
		wa->slices[i]=create_slice(wa->wave_data,start,nfft);
	}

	return wa;
}