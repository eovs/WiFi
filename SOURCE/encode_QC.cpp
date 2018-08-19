#include "stdafx.h"

#include <stdio.h>
#include <string.h>


static void unpackMatrix_double2int(double m[], int height, int width, int *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (int)m[j * height + i];
		}
	}
}

static void unpackRow_double2int(double m[], int height, int width, int result[] )
{
    int j;
	for( j = 0; j < width; ++j)
	{
		result[j] = (int)m[j];
	}
}

//% --------------------------------------------------%
//function y=cyclic_shift_left(x,s)
//% y is a cyclic shift of x by s positions left
//M=length(x);
//y=[x(s+1:M) x(1:s)];
static void cyclic_shift_left( int y[], int s, int yc[], int M )
{
   
    memcpy( yc, y + s, (M-s)*sizeof(y[0]) );
    memcpy( yc+(M-s), y, s*sizeof(y[0]) );
}



#define MAX_M 1000
void encode_qc( int codeword[], int *hb[], int M, int b, int c, int synd[] )
{
    int i, j;
    int k;
    int K = (c - b)*M;
	int L[3][2];
	int mode = 0;
	int r = M * b;
	static int u[MAX_M];
	static int y[MAX_M];
	static int yc[MAX_M];
	static int sumsynd[MAX_M];

	if( M > MAX_M )
		printf("WARNING M > MAX_M\n");

	j = 0;
	for( i = 0; i < b; i++ )
	{
		if( hb[i][c - b] >= 0 )
		{
			if( j > 2 )
				printf("weight of (c-b+1)st column is too large");
			else
			{
				L[j][0] = i;
				L[j][1] = hb[i][c - b];
				j++;
			}
		}
	}
	if( j < 3 )
		printf("weight of (c-b+1)st column is too small");

	if( L[0][1] == 0 && L[2][1] == 0 && L[1][1] > 0 )
		mode = 1; 
	
	if( L[0][1] != 0 && L[2][1] == L[0][1] && L[1][1] == 0 )
		mode = 2;

	if (mode == 0)
	{
		printf("(c-b+1)st column must be [0...A...0] or [A...0...A]");
		return;
	}

//%load message;
//codeword(1:k)=message;  // made before call
//% Compute partial syndrome
//synd=zeros(1,r);
//for j=1:c-b    % for all columns 
//    y=codeword((j-1)*M+(1:M)); % read block
//    for i=1:b  % for all rows
//        % circulate and add
//        if HB(i,j)>=0
//            yc=cyclic_shift_left(y,HB(i,j));
//            synd((i-1)*M+(1:M))=synd((i-1)*M+(1:M))+yc;
//        end;
//    end;
//end;
 
    memset( synd, 0, r * sizeof( synd[0] ) );
    for( j = 0; j < c-b; j++ )
    {
        // read block
        memcpy( y, codeword + j*M, M*sizeof( y[0] ) );  // M-items
        for( i = 0; i < b; i++ )                        // for all rows
        {
            //circulate and add
            if( hb[i][j] >= 0 )
            {
                cyclic_shift_left( y, hb[i][j], yc, M );    //yc=cyclic_shift_left(y,HB(i,j));
                //synd((i-1)*M+(1:M))=synd((i-1)*M+(1:M))+yc;
                for( k = 0; k < M; k++ )
                    synd[i*M +k] ^= yc[k];
            }
        }
    }
    //synd=mod(synd,2);     // done in cycle
    //% Compute sum of syndrom components
    //%sumsynd=mod(sum(reshape(synd,M,b),2),2)';
    //ss = reshape(synd,M,b);
    //sssum = sum( ss, 2 );
    //sumsynd = mod( sssum, 2 )';
    memset( sumsynd, 0, M*sizeof( sumsynd[0] ) );
    for( i = 0; i < b; i++ )
    {
        for( k = 0; k < M; k++ )
            sumsynd[k] ^= synd[i*M +k];
    }

	if( mode == 1 )
	{
		//% One check block is known
		//codeword(k+(1:M))=cyclic_shift_left(sumsynd,M-HB(L(2,1),c-b+1));
		cyclic_shift_left( sumsynd, M-hb[L[1][0]][c-b], codeword+K, M );
    
		//% Partial syndrom modification
		//synd((L(2,1)-1)*M+(1:M))=mod(synd((L(2,1)-1)*M+(1:M))+sumsynd,2);
		//synd((L(3,1)-1)*M+(1:M))=mod(synd((L(3,1)-1)*M+(1:M))+codeword(k+(1:M)),2);
		for( i = 0; i < M; i++ )
		{
			synd[(L[1][0])*M + i] ^= sumsynd[i];
			synd[(L[2][0])*M + i] ^= codeword[K+i];
		}

		//% recursion
		//for i=c-b+2:c
		//  j=i-c+b-2;
		//  codeword((i-1)*M+(1:M))=mod(synd(j*M+(1:M))+codeword((i-2)*M+(1:M)),2);
		//end;
		for( i = c-b+1; i < c; i++ )
		{
			j = i-c+b-1;
			for( k = 0; k < M; k++ )
				codeword[i*M+k] = synd[j*M+k] ^ codeword[(i-1)*M + k];
		}
	}

	if( mode == 2 )
	{
		//codeword(k+(1:M))=sumsynd;
		for( i = 0; i < M; i++ )
			codeword[K+i] = sumsynd[i];

		//u=cyclic_shift_left(sumsynd,HB(L(1,1),c-b+1)); % correction of syndromme 
		cyclic_shift_left( sumsynd, hb[L[0][0]][c-b], u, M );

		//synd((L(1,1)-1)*M+(1:M))=mod(synd((L(1,1)-1)*M+(1:M))+u,2);
		for( i = 0; i < M; i++ )
			synd[L[0][0]*M + i] ^= u[i];

		//synd((L(2,1)-1)*M+(1:M))=mod(synd((L(2,1)-1)*M+(1:M))+sumsynd,2);
		for( i = 0; i < M; i++ )
			synd[L[1][0]*M + i] ^= sumsynd[i];
				
		//synd((L(3,1)-1)*M+(1:M))=mod(synd((L(3,1)-1)*M+(1:M))+u,2);
		for( i = 0; i < M; i++ )
			synd[L[2][0]*M + i] ^= u[i];
	

		//% recursion
		//for i=c-b+2:c
		//	j=i-c+b-2;
		//	codeword((i-1)*M+(1:M))=synd(j*M+(1:M)); 
		//	synd((j+1)*M+(1:M))=mod(synd((j+1)*M+(1:M))+synd(j*M+(1:M)),2);
		//	%mod(synd(j*M+(1:M))+codeword((i-2)*M+(1:M)),2);
		//end;
		
		for( i = c-b+1; i < c; i++ )
		{
			j = i-c+b-1;
			for( k = 0; k < M; k++ )
			{
				codeword[i*M+k] = synd[j*M+k];
				synd[(j+1)*M+k] ^= synd[j*M+k];
			}
		}

	}


    //% CHECK for being a proper codeword
    //synd=zeros(1,r);
    //for j=1:c,
    //    y=codeword((j-1)*M+(1:M)); % read block
    //    for i=1:b,
    //        if HB(i,j)>=0
    //            yc=cyclic_shift_left(y,HB(i,j));
    //            synd((i-1)*M+(1:M))=synd((i-1)*M+(1:M))+yc;
    //        end;
    //    end; 
    //end;
    //synd=mod(synd,2);
    //if ~all(synd==0), error('bad coding'); end;
    memset( synd, 0, r*sizeof(synd[0]));
    for( j = 0; j < c; j++ )
    {
        memcpy( y, codeword + j*M, M*sizeof(y[0]));
        for( i = 0; i < b; i++ )
        {
            if( hb[i][j] >= 0 )
            {
                cyclic_shift_left( y, hb[i][j], yc, M );
                for( k = 0; k < M; k++ )
                {
                    synd[i*M +k] ^= yc[k];
                }
                 
            }
        }
    }
    //if ~all(synd==0), error('bad coding'); end;
    for( i = 0; i < r; i++ )
    {
        if( synd[i] )
//            mexErrMsgTxt("bad coding");
            printf("bad coding");
            
    }

}


