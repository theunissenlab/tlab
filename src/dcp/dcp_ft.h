#include <ss/ss.h>

#define MAX_SPEC_LEN 1024
#define MAX_SSPIKES 100    /* maximum number of sorted spikes */

typedef struct dcp_trial {
   char t_stamp[120];
	int nspikes;
	unsigned int *time;        /* time is in micro seconds */
	int npoints;
	short *waveform;
	SpikeList *ss;    
	int param;
	} DCP_TRIAL;

struct SEG_STRUCT {
   int begin;
   int end;         /* last element - location of cut */
   int len;         /* number of elements in the array */
   int subbegin;
   int subend;
   int sublen;     /* number of elements in each subsegment array */
};
typedef struct SEG_STRUCT seg_struct;


#define MAX_TONE_COMBS 50

typedef struct dcp_stim {
	char type[32];
	char songpath[128];
	char songfile[128];
	char segfile[128];
	char segcmd[128];
	int samprate;
	int dur;
	int length;
	int tone_freq;
	int sweep_f1;
	int sweep_f2;
	int nfreqs;
	int freq_array[MAX_TONE_COMBS];
	int rise_time;
	int fall_time;
	int nsegs;
   seg_struct *song_seg;
	short *stimulus;

  } DCP_STIM;



typedef struct dcp_data {
	char uv_version[120];
   char stim_spec[MAX_SPEC_LEN];
	char col_spec[MAX_SPEC_LEN];
	char prot_spec[MAX_SPEC_LEN];
	int n;
	int pre;
	int cd;
	int ici;
	int rand_ici;
	int adfreq;
	int adgain;
	int waveform;
	float prot_iscale;
	float prot_vscale;
	int ssflag;
	int prot_ra_num;
	float prot_ra_cmd;
   float prot_ra_start;
	float prot_ra_dur;
	float prot_ra_delay;
	int spike_list[MAX_SSPIKES]; 
	SpikeSet *spike_set;
	DCP_TRIAL *ptrial;
	DCP_STIM *pstim;
} DCP_DATA;

	
#define RET_ERR -1
#define RET_OK 0

DCP_STIM *setup_stimulus();

#define nint(x) ((int)(rint(x)))
