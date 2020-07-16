#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
char stim_spec[MAX_SPEC_LEN];
char file_name[256];
char cbuff[1024];
char temp_cbuff[1024];
DCP_DATA dcp_out;
DCP_TRIAL *ptout;
char ret_char;
char *str_item;
FILE *fin;
int i,j;

	if ( argc != 4 )
	{
		printf("Error: asci_to_dcp takes three arguments\n");
		printf("dcp_newpath <ascii_file> <sound_path> <sound_file>\n");
		exit(1);
	}

   /* Read file */
	fin = fopen(argv[1],"r");
	if ( fin == NULL ) 
	{
		printf("Error could not open input ascii file %s\n", argv[1]);
		exit(1);
	}

	/* Count lines in file */
	dcp_out.n = 0;
	while ( fgets(cbuff,1024,fin) ) (dcp_out.n)++;
	rewind(fin);

	/* Set dcp parameters */
   strcpy(dcp_out.uv_version,"ascii");
   sprintf(dcp_out.stim_spec,"stim song p=%s %s", argv[2], argv[3]);
	strcpy(dcp_out.col_spec,"unknown");
	dcp_out.pre = 0;
	dcp_out.cd = 3000;
	dcp_out.ici = 1000;
	dcp_out.rand_ici = 0;
	dcp_out.adfreq = 32000;
	dcp_out.adgain =1;
	dcp_out.waveform = 0;
	dcp_out.prot_iscale = 1.0;
	dcp_out.prot_vscale = 1.0;
	dcp_out.ssflag = 0;
	dcp_out.prot_ra_num = 0;
	dcp_out.prot_ra_cmd = 0.0;
	dcp_out.prot_ra_start = 0.0;
	dcp_out.prot_ra_dur = 0.0;
	dcp_out.prot_ra_delay = 0.0;
	dcp_out.spike_set = NULL; 
	dcp_out.pstim  = NULL; 
	dcp_out.ptrial = (DCP_TRIAL *)calloc(dcp_out.n, sizeof(DCP_TRIAL));


	for ( i = -1; ++i < dcp_out.n; )
	{
		ptout = dcp_out.ptrial + i;
		fgets(cbuff,1024,fin);
		cbuff[strlen(cbuff)-1]=0; 
		printf("%s\n",cbuff);
		strcpy(temp_cbuff,cbuff);

		str_item = strtok(cbuff," \t");
		j = 0;
		while ( str_item != NULL )
		{
			j++;
			str_item = strtok(NULL," \t");
		}
		ptout->nspikes=j;
		ptout->time = calloc(j,sizeof(int));
		ptout->npoints=0;
		ptout->waveform=NULL;
		ptout->ss=NULL;
		ptout->param=NULL;

		str_item = strtok(temp_cbuff," \t");
		j = 0;
		printf("Trial %d (%d spikes):", i, ptout->nspikes);
		while ( str_item != NULL )
		{
			ptout->time[j++]=(int)(atof(str_item));
			str_item = strtok(NULL," \t");
			printf("%d ",ptout->time[j-1]);
		}
		printf("\n");
	}

	sprintf(file_name,"%s.dcp",argv[1]);			
   if ( write_dcp_data(&dcp_out, file_name) == RET_ERR )
		printf("Error writing out dcp file\n");

}

