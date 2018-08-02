#include "stdafx.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "modulation.h"


static int Qav[4] = { 4, 16, 64, 256 };

int findQ( int q )
{
    for( int i = 0; i < 4; i++ )
    {
        if( q == Qav[i] )
            return i;
    }
    return NOT_FOUND;
                
}


QAM_MODULATOR_STATE* QAM_modulator_open( int Q, int L, int m )
{
    QAM_MODULATOR_STATE* st;
    int i,j;
    
    st = (QAM_MODULATOR_STATE*)calloc(sizeof(QAM_MODULATOR_STATE),1);
    if( !st ) return NULL;
    
    st->Lorg = L;
    st->Q = Q;
    st->m = m;
    
    int r = L % m;
    if( r )
        st->Lfact = L+m-r;
    else
        st->Lfact = L;
    st->ns = st->Lfact / m;
    
    st->x = (short*)calloc( st->Lfact, sizeof(st->x[0]) );
    if( !st->x ) return NULL;
    
    st->z1 = (int*)calloc( st->ns, sizeof(st->z1[0]));
    if( !st->z1 ) return NULL;
    st->z2 = (int*)calloc( st->ns, sizeof(st->z2[0]));
    if( !st->z2 ) return NULL;
    
    for( i = m/2-1, j = 0; i>= 0; i--, j++ )
        st->p[j] = 1<<i;
    
    st->dx = (int*)calloc( st->Lfact, sizeof(st->dx[0]) );
    return st;
    
}

void QAM_modulator_close(QAM_MODULATOR_STATE* st)
{
    st->Lorg = 0;
    st->Q = 0;
    st->Lfact = 0;
    free( st->x );
    free( st->z1 );
    free( st->z2 );
	free( st->dx );
    
    free( st );
}


static short gray[]={0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13,  8,  9, 11, 10};  //% anti-gray
static short s[] = {0,1,3,7, 15};

static void GrayPAM( int x[], int size, int y[], int order )
{
 
    for( int i = 0; i < size; i++ )
         y[i] = 2*gray[x[i]]-s[order];
 }


void QAM_modulator( QAM_MODULATOR_STATE *qam_mod_state, int *in, int *out[2]  )
{
	int n, i, j, l;
	int m = qam_mod_state->m;
    int mdiv2 = qam_mod_state->m/2;
   
	for( n = 0, j = 0; n < qam_mod_state->Lfact; n += qam_mod_state->m, j++ )
    {
        int sum = 0;
        for( i = 0; i < mdiv2; i++ )
            sum += (int)(qam_mod_state->p[i] * in[n+i]);

        qam_mod_state->z1[j] = sum;
        
        sum = 0;
        for( i = mdiv2, l = 0; i < m; i++, l++ )
            sum += (int)(qam_mod_state->p[l] * in[n+i]);
 
        qam_mod_state->z2[j] = sum;
        
    }

    GrayPAM( qam_mod_state->z1, qam_mod_state->ns, out[0], mdiv2 );
    GrayPAM( qam_mod_state->z2, qam_mod_state->ns, out[1], mdiv2 );
}



