#ifndef _PARAMS_H_
#define _PARAMS_H_

typedef struct
{
	int cfg;
	int n_pass;
	int n_iter;
	int n_event;
	int llr_bits;
	int code_M;
	int dec_type;
	double target_FER;
	double SNR_start;
	double SNR_stop;
	double SNR_step;
	int SNR_flag;
	int llr_limit;
	char file_name[256];
	int n_attempts;
} SIMULATION_PARAMS;

int set_params( char *fileName, SIMULATION_PARAMS *params );

#endif //_PARAMS_H_
