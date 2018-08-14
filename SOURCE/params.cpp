#include "stdafx.h"

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
	{ "target_FER", 7 }
};

int set_params( char *fileName, SIMULATION_PARAMS *params )
{
	char line[1000];
	char name[100];
	char val[100];
	int nparams = sizeof(supported_params) / sizeof(supported_params[0]);
	int i;

	FILE *fp;

	fopen_s( &fp, fileName, "rt" );
	if( fp == NULL )
		return 0;

	printf("CONFIG FILE:\n");
	while( !feof( fp ) )
	{
		char *ptr = fgets( line, sizeof(line), fp );
		if( ptr == NULL )
			break;

		sscanf_s( line, "%s%s", name, (int)sizeof(name), val, (int)sizeof(val) );

		for( i = 0; i < nparams; i++ )
		{
			if( strcmp( name, supported_params[i].paramName ) == 0 )
			{
				printf("%s", line);
		
				switch( supported_params[i].paramIndex )
				{
				case 0: params->cfg        = atoi( val );	break;
				case 1: params->n_pass     = atoi( val );	break;
				case 2: params->n_iter     = atoi( val );	break;
				case 3: params->n_event    = atoi( val );	break;
				case 4: params->llr_bits   = atoi( val );	break;
				case 5: params->code_M     = atoi( val );	break;
				case 6: params->dec_type   = atoi( val );	break;
				case 7: params->target_FER = atof( val );	break;
				}

				break;
			}
		}

		//if( i == nparams )	printf("parametr %s unsupported\n", name );

	}

	fclose( fp );
	return 1;
}
