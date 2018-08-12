#include "stdafx.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "modulation.h"


static int Qav[4] = { 4, 16, 64, 256 };

static int findQ( int q )
{
    for( int i = 0; i < 4; i++ )
    {
        if( q == Qav[i] )
            return i;
    }
    return NOT_FOUND;
                
}

QAM_DEMODULATOR_STATE* QAM_demodulator_open( double T, double sigma, int Q, int n, int m, int ns, int out_type )
{
    QAM_DEMODULATOR_STATE* st;
    
    st = (QAM_DEMODULATOR_STATE*)calloc(1, sizeof(QAM_DEMODULATOR_STATE));
    if( !st ) return NULL;
    
    st->T = T;
    st->sigma = sigma;
    st->Q = Q;
    st->n = n;
    st->m = m;
    st->ns = ns;
    st->DemodOutType = out_type;
  
    return st;
}

void QAM_demodulator_close(QAM_DEMODULATOR_STATE* st )
{
    free( st );
}


void Demodulate( QAM_DEMODULATOR_STATE* st, double *pMod[2], double pRes[] )
{
    int m = st->m;
    int i;
    int ns = st->ns;
    int n = st->n;
    double sigma = st->sigma;
    double T = st->T;
    
   
    if( m == 2 )
    {   
        int i,j;
        double P;
        double sigma2 = sigma*sigma;
        for( i = j = 0; i < n; i+=2, j++ )
            pRes[i] = 2.0*pMod[0][j]/sigma2;
        for( i = 1, j = 0; i < n; i+=2, j++ )
            pRes[i] = 2.0*pMod[1][j]/sigma2;
        if( st->DemodOutType )
        {
            P= 0.0;
            for( i = 0; i < n; i++ )
            {
                pRes[i] = exp( pRes[i] );
                P += pRes[i];
            }
            for( i = 0; i< n; i++ )
                pRes[i] /= P;
        }
       return;
    }

//N0=2*sigma^2;
//SQ=2^(m/2);   % square root of Q; 
    double N0 = 2.0 * sigma *sigma;
    int SQ = 1 <<( m/2 );   //square root of Q;


    //s=[1,3,7,15];
    //lattice=2*(0:(SQ-1))-s(m/2);
    static int s[] = { 0,1,3,7,15 };
    static int lattice[16];
    
    for( int i = 0; i < SQ; i++ )
    {
        lattice[i] = 2*i - s[m/2];
    }
//    mexPrintf("\n");

//L=zeros(m,ns);
//P1=zeros(m,ns);
//P=zeros(1,SQ);
//    memset( st->L, 0, n*sizeof(st->L[0]));
//    memset( st->P1, 0, n*sizeof(st->P1[0]));
    double P[16] = { 0 };
//    double D[16];
    
//    for( int i = 0; i < ns; i++ )
    int  ix;
    for( i = ix = 0; i < n; i+=m, ix++ )
    {
        int h = -1; //0;
        for( int j = 0; j < 2; j++ )
        {
             //D=(x(j,i)-lattice).^2/N0;
             //P(D<T)=exp(-D(D<T));
             //P(D>=T)=0;
            double sum = 0;
            for( int i1 = 0; i1 < SQ; i1++ )
            {
                double tmp = pMod[j][ix] - lattice[i1];
                tmp *= tmp;
                tmp /= N0;
  //              D[i1] = tmp;
                if( tmp < T )
                    P[i1] = exp(-tmp);
                else
                    P[i1] = 0.0;
                sum +=  P[i1];
            }

            //P=P/sum(P);
            //double sum = 0;
            //for( int i1 = 0; i1 < SQ; i1++)
            //    sum += P[i1];
            for( int i1 = 0; i1 < SQ; i1++ )
            {
                P[i1] /= sum;
            }

            switch( m/2 )
            {
            case 2: //% QAM-16
            {  // %  00 01 11 10 
                //   p0=P(1)+P(2);
                //   p1=P(3)+P(4);
                //   h=h+1;
                double p0 = P[0] + P[1];
                double p1 = P[2] + P[3];
                h++;
               //if p0==0, L(h,i)=T; P1(h,i)=1; 
               //else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
               //     else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
               //end; 
                if( p0 == 0.0 )
                {
                    
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1 == 0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = -T;
                        else
                            pRes[h/*ns*/+i] = 0.0;
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1;
                    }
                }
                //   p0=P(1)+P(4);
                //   p1=P(2)+P(3);
                //   h=h+1;
                p0 = P[0] + P[3];
                p1 = P[1] + P[2];
                h++;
                //if p0==0, L(h,i)=T; P1(h,i)=1; 
                //else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
                //     else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
                //end;      
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1 == 0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = -T;
                        else
                            pRes[h/*ns*/+i] = 0.0;
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1;
                    }
                }
                break;
            }
            case 3: // % QAM-64 
            {
           
                //p12=P(1)+P(2);
                //p34=P(3)+P(4);
                //p56=P(5)+P(6);
                //p78=P(7)+P(8);
                //p1234=p12+p34;
                //p5678=p56+p78;
                //p1278=p12+p78;
                //p3456=p34+p56;
                //p0=p1234; p1=p5678;
                //=h+1;
                double p12 = P[0]+P[1];
                double p34 = P[2]+P[3];
                double p56 = P[4]+P[5];
                double p78 = P[6]+P[7];
                double p1234 = p12+p34;
                double p5678 = p56+p78;
                double p1278 = p12+p78;
                double p3456 = p34+p56;
                double p0 = p1234;
                double p1 = p5678;
                h++;
//                if p0==0, L(h,i)=T; P1(h,i)=1; 
//                 else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//                     else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//                 end;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
//           p0=p1278; p1=p3456;
//           h=h+1;
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                p0 = p1278; p1 = p3456;
                h++;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
//           p0=P(1)+P(4)+P(5)+P(8);
//           p1=P(2)+P(3)+P(6)+P(7);
//           h=h+1;
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                p0 = P[0] + P[3] + P[4] + P[7];
                p1 = P[1] + P[2] + P[5] + P[6];
                h++;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
            } // case 3
            break;
            case  4:
            {
//           p12=P(1)+P(2);
//           p34=P(3)+P(4);
//           p56=P(5)+P(6);
//           p78=P(7)+P(8);
//           p9A=P(9)+P(10);
//           pBC=P(11)+P(12);
//           pDE=P(13)+P(14);
//           pFG=P(15)+P(16);
//           p1234=p12+p34;
//           p5678=p56+p78;
//           p9ABC=p9A+pBC;
//           pDEFG=pDE+pFG;
//           p1to8=p1234+p5678;
//           p9toG=p9ABC+pDEFG;
//           p0=p1to8; p1=p9toG;
//           h=h+1;
                double p12 = P[0] + P[1];
                double p34 = P[2] + P[3];
                double p56 = P[4] + P[5];
                double p78 = P[6] + P[7];
                double p9A = P[8] + P[9];
                double pBC = P[10]+P[11];
                double pDE = P[12]+P[13];
                double pFG = P[14]+P[15];
                double p1234 = p12 + p34;
                double p5678 = p56 + p78;
                double p9ABC = p9A + pBC;
                double pDEFG = pDE + pFG;
                double p1to8 = p1234 + p5678;
                double p9toG = p9ABC + pDEFG;
                double p0 = p1to8;
                double p1 = p9toG;
                h++;    //0
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
//           p1=p5678+p9ABC;
//           p0=p1234+pDEFG;
//           h=h+1;
                p1 = p5678 + p9ABC;
                p0 = p1234 + pDEFG;
                h++;    //1
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
//           p1=p34+p56+pBC+pDE;
//           p0=p12+p78+p9A+pFG;
//           h=h+1;
                p1 = p34 + p56 + pBC + pDE;
                p0 = p12 + p78 + p9A + pFG;
                h++;    //2
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
//           p1=P(2)+P(3)+P(6)+P(7)+P(10)+P(11)+P(14)+P(15);
//           p0=P(1)+P(4)+P(5)+P(8)+P( 9)+P(12)+P(13)+P(16);
//           h=h+1;
                p1 = P[1]+P[2]+P[5]+P[6]+P[9]+P[10]+P[13]+P[14];
                p0 = P[0]+P[3]+P[4]+P[7]+P[8]+P[11]+P[12]+P[15];
                h++;    //3
//           if p0==0, L(h,i)=T; P1(h,i)=1; 
//           else if p1==0,L(h,i)=-T;  P1(h,i)=0; 
//               else  L(h,i)=log(p1/p0); P1(h,i)=p1; end; 
//           end;
                if( p0 == 0.0 )
                {
                    if( st->DemodOutType == 0 )
                        pRes[h/*ns*/+i] = T;
                    else
                        pRes[h/*ns*/+i] = 1.0;
                }
                else
                {
                    if( p1==0.0 )
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i]=-T;
                        else
                            pRes[h/*ns*/+i] = 0.0; 
                    }
                    else
                    {
                        if( st->DemodOutType == 0 )
                            pRes[h/*ns*/+i] = log(p1/p0);
                        else
                            pRes[h/*ns*/+i] = p1; 
                    }
                }
            } // case 4
            break;
            }// switch
        }
    }
}
