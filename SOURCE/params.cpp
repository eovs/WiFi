#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"

typedef struct  
{
	char paramName[100];
	int paramIndex;
} SUPPORTED_PARAMS;

SUPPORTED_PARAMS supported_params[] = 
{
	{ "cfg",        0 },        
	{ "n_pass",     1 },    
	{ "n_iter",     2 },    
	{ "n_event",    3 },   
	{ "llr_bits",   4 },   
	{ "code_M",     5 },    
	{ "dec_type",   6 },  
	{ "target_FER", 7 },
	{ "snr",        8 },
	{ "aux_file",   9 }
};

char* get_string( char line[], char *val )
{
	char string[100];
	size_t len = 0;

	if( sscanf_s( line, "%s", string, (int)sizeof(string) ) == 1)
	{
		strcpy_s( val,  (int)sizeof(string), string );
		len = strlen( string );
	}
	return  strstr( line, string ) + len;
}

char* get_int( char line[], int *val )
{
	char string[100];
	size_t len = 0;

	if( sscanf_s( line, "%s", string, (int)sizeof(string) ) == 1)
	{
		*val = atoi( string );
		len = strlen( string );
	}
	return  strstr( line, string ) + len;
}

char* get_dbl( char line[], double *val )
{
	char string[100];
	size_t len = 0;

	if( sscanf_s( line, "%s", string, (int)sizeof(string) ) == 1)
	{
		*val = atof( string );
		len = strlen( string );
	}
	return  strstr( line, string ) + len;
}

int set_params( char *fileName, SIMULATION_PARAMS *params )
{
	char line[1000];
	char name[100];
	int nparams = sizeof(supported_params) / sizeof(supported_params[0]);
	int i;

	FILE *fp;

	fopen_s( &fp, fileName, "rt" );
	if( fp == NULL )
		return 0;

	printf("CONFIG FILE:\n");
	params->SNR_flag = 0;

	while( !feof( fp ) )
	{
		char *curr_line = fgets( line, sizeof(line), fp );
		if( curr_line == NULL )
			break;
		
		if( line[0] == '\n' )
			break;

		curr_line = strstr( line, "//" ); 
		curr_line[0] = '\n';
		curr_line[1] = '\0';

		curr_line = get_string( line, name );
		if( curr_line == NULL ) continue;
		
		for( i = 0; i < nparams; i++ )
		{
			if( strcmp( name, supported_params[i].paramName ) == 0 )
			{
				printf("%s", line);
		
				switch( supported_params[i].paramIndex )
				{
				case 0: curr_line = get_int( curr_line, &params->cfg );			break;
				case 1: curr_line = get_int( curr_line, &params->n_pass );		break;
				case 2: curr_line = get_int( curr_line, &params->n_iter );		break;
				case 3: curr_line = get_int( curr_line, &params->n_event );		break;
				case 4: curr_line = get_int( curr_line, &params->llr_bits );	break;
				case 5: curr_line = get_int( curr_line, &params->code_M );		break;
				case 6: curr_line = get_int( curr_line, &params->dec_type );	break;
				case 7: curr_line = get_dbl( curr_line, &params->target_FER );	break;
				case 8: curr_line = get_dbl( curr_line, &params->SNR_start );
					    curr_line = get_dbl( curr_line, &params->SNR_stop );
						curr_line = get_dbl( curr_line, &params->SNR_step );	
						params->SNR_flag  = 1; break;
				case 9: curr_line = get_string( curr_line, params->file_name );	break;
				}

				break;
			}
		}

		//if( i == nparams )	printf("parametr %s unsupported\n", name );

	}

	fclose( fp );
	return 1;
}
