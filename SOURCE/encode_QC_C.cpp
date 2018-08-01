#include "stdafx.h"

#ifndef SKIP_MEX
#include <mex.h>
#endif

#include <string.h>


void unpackMatrix_double2int(double m[], int height, int width, int *result[] )
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

void unpackRow_double2int(double m[], int height, int width, int result[] )
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
void cyclic_shift_left( int y[], int s, int yc[], int M )
{
   
    memcpy( yc, y + s, (M-s)*sizeof(y[0]) );
    memcpy( yc+(M-s), y, s*sizeof(y[0]) );
}

#ifndef SKIP_MEX
int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;

	p = (int**)calloc( b, sizeof(int*) );
    if( !p )    
        mexErrMsgTxt("Allocation error");
    p[0] = (int*)calloc(b*c, sizeof(int));
    if( !p[0] )
        mexErrMsgTxt("Allocation error");
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c; //(int*)calloc( c, sizeof(int) );
	}
	return p;
}
void free2d_int( int **p )
{
	int i;

//	for( i = 0; i < b; i++ )
//		free( p[i] );
    free( p[0] );
	free( p );
}
#endif

#define MAX_M 1000
void encode_qc( int codeword[], short *hb[], int M, int b, int c, short synd[] )
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


#ifndef SKIP_MEX
// %
// % encode_QC_C interface:
// % 1. SETUP ( before simulation )
// % ========
// % result = encode_QC_C( HB, M );
// %  result == 1, if OK
// %
// % 2. DECODE
// % =========
// %  codeword = encode_QC_C( message );
// %     codeword  -  enccoding result
// %
// % 3. CLEAR ( free memory allocations )
// % ========
// %  encode_QC_C();     
// % 

void mexFunction(int nOut, mxArray *pOut[], int nInp, const mxArray *pInp[])
{
    static int **hb = NULL;         // b x c
    static int *synd = NULL;        // r
    static int *codeword = NULL;    // n
    static int *y;                  // M
    static int *yc;                 // M
    static int *sumsynd;            // M
    static int b, c, r, n, k, M;
    static int L[3][2] = { 0 };
    int i, j;
    
    if( nInp == 2 )
    {
        // SETUP :result = encode_QC_C( HB, M );
        if( nOut != 1 )
        {
            mexErrMsgTxt("Only one output argument allowed for ENCODING SETUP");
        }

        // Setup call for HB convertion to C-style and for store in dynamic mamory
        // Correct call:  encode_QC_C( HB, M);
        
        M = (int)mxGetPr(pInp[1])[0];
        
        //[b,c]=size(HB);

        b = (int)mxGetM(pInp[0]);
        c = (int)mxGetN(pInp[0]);

        r=b*M;   //% number of rows
        n=c*M;   //% number of columns
        k=n-r;   //% message length  
        
        // Check if dynamic arrays were allocated
        if( hb ) free2d_int( hb );
        if( codeword ) free(codeword);
        if( synd ) free( synd );
        if( y ) free( y );
        if( yc ) free( yc );
        if( sumsynd ) free( sumsynd );
        
        
        hb = Alloc2d_int( b, c );
        if( !hb ) mexErrMsgTxt("Allocation error on ENCODING SETUP");

        codeword = (int*)calloc( n, sizeof(codeword[0]) );
        if( !codeword ) mexErrMsgTxt("Allocation error on ENCODING SETUP");
        
        synd = (int*)calloc( r, sizeof(synd[0]) );
        if( !synd ) mexErrMsgTxt("Allocation error on ENCODING SETUP");

        y = (int*)calloc( M, sizeof( y[0] ) );
        if( !y ) mexErrMsgTxt("Allocation error on ENCODING SETUP");

        yc = (int*)calloc( M, sizeof( yc[0] ) );
        if( !yc ) mexErrMsgTxt("Allocation error on ENCODING SETUP");
        
        sumsynd = (int*)calloc( M, sizeof( sumsynd[0] ) );
        if( !sumsynd ) mexErrMsgTxt("Allocation error on ENCODING SETUP");
        
        
        unpackMatrix_double2int(mxGetPr(pInp[0]), b, c, hb );

#if 0        
        //j=0;
        //for i=1:b 
        //  if HB(i,c-b+1)>=0,
        //      j=j+1;
        //      if j>3, error('weight of (c-b+1)st column is too large'); end;
        //      L(j,1)=i;
        //      L(j,2)=HB(i,c-b+1);
        //  end;
        //end;
        //if j<3, error('weight of (c-b+1)st column is too small'); end;
        //if L(1,2)~=0 || L(3,2)~=0 || L(2,2)==0,  
        //  error('(c-b+1)st column must be [0...A...0]');
        //end;
#endif   
        j = 0;
        for( i = 0; i < b; i++ )
        {
            if( hb[i][c-b] >= 0 )
            {
                j++;
                if( j > 3 )
                    mexErrMsgTxt("weight of (c-b+1)st column is too large");
                L[j-1][0] = i+1;
                L[j-1][1] = hb[i][c-b];
            }
        }
        if( j < 3 )
            mexErrMsgTxt("weight of (c-b+1)st column is too small");
        if( L[0][1] != 0 || L[2][1] != 0 || L[1][1] == 0 )
            mexErrMsgTxt("(c-b+1)st column must be [0...A...0]");
        //mexPrintf("nh = %d mh = %d N = %d maxiter = %d\n", nh, mh, N, maxiter);
        //mexPrintf("%d %d\n", L[0][0], L[0][1] );
        //mexPrintf("%d %d\n", L[1][0], L[1][1] );
        //mexPrintf("%d %d\n", L[2][0], L[2][1] );
        
        pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
        *(mxGetPr(pOut[0])) = 1;
        return; // end of setup
    }
    
    #if 1
    if( nInp == 1 )
    {
        double *p;
        if( nOut != 1 )
            mexErrMsgTxt("Only one output argument allowed for ENCODING");
        // Decoding, only 1 input param - input vector
        unpackRow_double2int(mxGetPr(pInp[0]), 1, k, codeword );
        memset( codeword + k, 0, r*sizeof( codeword[0] ) );
        //iter = bp_decod_qc_lm( y, hard, hd, nh, mh, M, maxiter);
        encode_qc( codeword, hb, M, b, c, r, synd, y, yc, sumsynd, L );

        pOut[0]=mxCreateDoubleMatrix(1,n,mxREAL);
        p = mxGetPr(pOut[0]);
        for( i = 0; i < n; i++ )
            p[i] = codeword[i];
        //pOut[1] =  mxCreateDoubleScalar(iter);
        //pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
        //*(mxGetPr(pOut[0])) = iter;
        return;
    }
    #endif
    
    if( nInp == 0 )
    {
        // Check if dynamic arrays were allocated
        if( hb )
            free2d_int( hb );
        hb = NULL;
        if( codeword )
            free(codeword);
        codeword = NULL;
        if( synd )
            free(synd);
        synd = NULL;
        if( y )
            free( y );
        y = NULL;
        if( yc )
            free( yc );
        yc = NULL;
        if( sumsynd )
            free( sumsynd );
        sumsynd = NULL;
    }
}
#endif