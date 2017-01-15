/***************************************************************
**                                                            **
**        D I F F E R E N T I A L     E V O L U T I O N       **
**                                                            **
** Program: de.c                                              **
** Version: 3.6                                               **
**                                                            **
** Authors: Dr. Rainer Storn                                  **
**          c/o ICSI, 1947 Center Street, Suite 600           **
**          Berkeley, CA 94707                                **
**          Tel.:   510-642-4274 (extension 192)              **
**          Fax.:   510-643-7684                              **
**          E-mail: storn@icsi.berkeley.edu                   **
**          WWW: http://http.icsi.berkeley.edu/~storn/        **
**          on leave from                                     **
**          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6          **
**          D-81739 Muenchen, Germany                         **
**          Tel:    636-40502                                 **
**          Fax:    636-44577                                 **
**          E-mail: rainer.storn@zfe.siemens.de               **
**                                                            **
**          Kenneth Price                                     **
**          836 Owl Circle                                    **
**          Vacaville, CA 95687                               **
**          E-mail: kprice@solano.community.net               ** 
**                                                            **
** This program implements some variants of Differential      **
** Evolution (DE) as described in part in the techreport      **
** tr-95-012.ps of ICSI. You can get this report either via   **
** ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z  **
** or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
** A more extended version of tr-95-012.ps is submitted for   **
** publication in the Journal Evolutionary Computation.       ** 
**                                                            **
** You may use this program for any purpose, give it to any   **
** person or change it according to your needs as long as you **
** are referring to Rainer Storn and Ken Price as the origi-  **
** nators of the the DE idea.                                 **
** If you have questions concerning DE feel free to contact   **
** us. We also will be happy to know about your experiences   **
** with DE and your suggestions of improvement.               **
**                                                            **
***************************************************************/
/**H*O*C**************************************************************
**                                                                  **
** No.!Version! Date ! Request !    Modification           ! Author **
** ---+-------+------+---------+---------------------------+------- **
**  1 + 3.1  +5/18/95+   -     + strategy DE/rand-to-best/1+  Storn **
**    +      +       +         + included                  +        **
**  1 + 3.2  +6/06/95+C.Fleiner+ change loops into memcpy  +  Storn **
**  2 + 3.2  +6/06/95+   -     + update comments           +  Storn **
**  1 + 3.3  +6/15/95+ K.Price + strategy DE/best/2 incl.  +  Storn **
**  2 + 3.3  +6/16/95+   -     + comments and beautifying  +  Storn **
**  3 + 3.3  +7/13/95+   -     + upper and lower bound for +  Storn **
**    +      +       +         + initialization            +        **
**  1 + 3.4  +2/12/96+   -     + increased printout prec.  +  Storn **
**  1 + 3.5  +5/28/96+   -     + strategies revisited      +  Storn **
**  2 + 3.5  +5/28/96+   -     + strategy DE/rand/2 incl.  +  Storn **
**  1 + 3.6  +8/06/96+ K.Price + Binomial Crossover added  +  Storn **
**  2 + 3.6  +9/30/96+ K.Price + cost variance output      +  Storn **
**  3 + 3.6  +9/30/96+   -     + alternative to ASSIGND    +  Storn **
**  4 + 3.6  +10/1/96+   -    + variable checking inserted +  Storn **
**  5 + 3.6  +10/1/96+   -     + strategy indic. improved  +  Storn **
**                                                                  **
***H*O*C*E***********************************************************/
 
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "memory.h"
#include <time.h>
 
// 最大族群數, NP
#define MAXPOP  5000
// 最大向量維度, D
#define MAXDIM  35
#define MAXIMAPROBLEM 0
#define PENALITY 1000
 
/*------Constants for rnd_uni()--------------------------------------------*/
 
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
 
// 與機構合成相關的常數定義
#define PI 3.1415926
#define degree PI/180.0
#define mech_loop -1
#define NUM_OF_POINTS 10
 
/*------------------------Macros----------------------------------------*/
 
/*#define ASSIGND(a,b) memcpy((a),(b),sizeof(double)*D) */  /* quick copy by Claudio */
                                                           /* works only for small  */
                                                           /* arrays, but is faster.*/
 
/*------------------------Globals---------------------------------------*/
 
long  rnd_uni_init;                 /* serves as a seed for rnd_uni()   */
double c[MAXPOP][MAXDIM], d[MAXPOP][MAXDIM];
double (*pold)[MAXPOP][MAXDIM], (*pnew)[MAXPOP][MAXDIM], (*pswap)[MAXPOP][MAXDIM];
 
/*---------Function declarations----------------------------------------*/
 
void  assignd(int D, double a[], double b[]);
double rnd_uni(long *idum);    /* uniform pseudo random number generator */
double extern evaluate(int D, double tmp[], long *nfeval); /* obj. funct. */
 
// 與機構合成相關的函式宣告
double distance(double x0, double y0, double x1, double y1);
double rr(double L1, double dd, double theta);
struct Coord triangletip_coord( double x0, double y0, double R0, double R1, double x1, double y1, double localt);
void mechanism(double x0, double y0, double x1, double y1, double L1,
  double L2, double L3, double L5, double L6, double input_angles[NUM_OF_POINTS], struct Coord output_points[NUM_OF_POINTS]);
double error_function(struct Coord output_points[NUM_OF_POINTS], struct Coord target_points[NUM_OF_POINTS]);
   
struct Coord finaltip_coord(struct Coord tip1_coord, struct Coord tip2_coord, double r1, double r2);
 
/*---------Function definitions-----------------------------------------*/
// 指定向量 b 為 a
void  assignd(int D, double a[], double b[])
/**C*F****************************************************************
**                                                                  **
** Assigns D-dimensional vector b to vector a.                      **
** You might encounter problems with the macro ASSIGND on some      **
** machines. If yes, better use this function although it's slower. **
**                                                                  **
***C*F*E*************************************************************/
{
   int j;
   for (j=0; j<D; j++)
   {
      a[j] = b[j];
   }
}
 
// 產生 0 ~ 1 間的亂數
double rnd_uni(long *idum)
/**C*F****************************************************************
**                                                                  **
** SRC-FUNCTION   :rnd_uni()                                        **
** LONG_NAME      :random_uniform                                   **
** AUTHOR         :(see below)                                      **
**                                                                  **
** DESCRIPTION    :rnd_uni() generates an equally distributed ran-  **
**                 dom number in the interval [0,1]. For further    **
**                 reference see Press, W.H. et alii, Numerical     **
**                 Recipes in C, Cambridge University Press, 1992.  **
**                                                                  **
** FUNCTIONS      :none                                             **
**                                                                  **
** GLOBALS        :none                                             **
**                                                                  **
** PARAMETERS     :*idum    serves as a seed value                  **
**                                                                  **
** PRECONDITIONS  :*idum must be negative on the first call.        **
**                                                                  **
** POSTCONDITIONS :*idum will be changed                            **
**                                                                  **
***C*F*E*************************************************************/
{
  long j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
 
  if (*idum <= 0)
  {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--)
    {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
 
}/*------End of rnd_uni()--------------------------*/
 
// 將上下限轉為全域變數
double inibound_h;      /* upper parameter bound              */
double inibound_l;      /* lower parameter bound              */
// 與機構合成相關的全域變數
// 宣告一個座標結構
struct Coord {
    double x;
    double y;
  // 這裡保留 double z;
};
 
main(int argc, char *argv[])
/**C*F****************************************************************
**                                                                  **
** SRC-FUNCTION   :main()                                           **
** LONG_NAME      :main program                                     **
** AUTHOR         :Rainer Storn, Kenneth Price                      **
**                                                                  **
** DESCRIPTION    :driver program for differential evolution.       **
**                                                                  **
** FUNCTIONS      :rnd_uni(), evaluate(), printf(), fprintf(),      **
**                 fopen(), fclose(), fscanf().                     **
**                                                                  **
** GLOBALS        :rnd_uni_init    input variable for rnd_uni()     **
**                                                                  **
** PARAMETERS     :argc            #arguments = 3                   **
**                 argv            pointer to argument strings      **
**                                                                  **
** PRECONDITIONS  :main must be called with three parameters        **
**                 e.g. like de1 <input-file> <output-file>, if     **
**                 the executable file is called de1.               **
**                 The input file must contain valid inputs accor-  **
**                 ding to the fscanf() section of main().          **
**                                                                  **
** POSTCONDITIONS :main() produces consecutive console outputs and  **
**                 writes the final results in an output file if    **
**                 the program terminates without an error.         **
**                                                                  **
***C*F*E*************************************************************/
 
{
   char  chr;             /* y/n choice variable                */
   char  *strat[] =       /* strategy-indicator                 */
   {
            "",
            "DE/best/1/exp",
            "DE/rand/1/exp",
            "DE/rand-to-best/1/exp",
            "DE/best/2/exp",
            "DE/rand/2/exp",
            "DE/best/1/bin",
            "DE/rand/1/bin",
            "DE/rand-to-best/1/bin",
            "DE/best/2/bin",
            "DE/rand/2/bin"
   };
 
   int   i, j, L, n;      /* counting variables                 */
   int   r1, r2, r3, r4;  /* placeholders for random indexes    */
   int   r5;              /* placeholders for random indexes    */
   int   D;               /* Dimension of parameter vector      */
   int   NP;              /* number of population members       */
   int   imin;            /* index to member with lowest energy */
   int   refresh;         /* refresh rate of screen output      */
   int   strategy;        /* choice parameter for screen output */
   int   gen, genmax, seed;   
 
   long  nfeval;          /* number of function evaluations     */
 
   double trial_cost;      /* buffer variable                    */
   // 將上下限轉為全域變數, 可能要根據各變數加以設定
   //double inibound_h;      /* upper parameter bound              */
   //double inibound_l;      /* lower parameter bound              */
   double tmp[MAXDIM], best[MAXDIM], bestit[MAXDIM]; /* members  */
   double cost[MAXPOP];    /* obj. funct. values                 */
   double cvar;            /* computes the cost variance         */
   double cmean;           /* mean cost                          */
   double F,CR;            /* control variables of DE            */
   double cmin;            /* help variables                     */
 
   FILE  *fpin_ptr;
   FILE  *fpout_ptr;
    
// 計算執行過程所需時間起點, 需要導入 time.h
  clock_t start = clock();
   
/*------Initializations----------------------------*/
 
 //if (argc != 3)                                 /* number of arguments */
 //{
    //printf("\nUsage : de <input-file> <output-file>\n");
    //exit(1);
 //}
 
// 將結果寫入 out.dat
 fpout_ptr = fopen("out.dat","w");          /* open output file for reading,    */
                                          /* to see whether it already exists */
 /*
 if ( fpout_ptr != NULL )
 {
    printf("\nOutput file %s does already exist, \ntype y if you ",argv[2]);
    printf("want to overwrite it, \nanything else if you want to exit.\n");
    chr = (char)getchar();
    if ((chr != 'y') && (chr != 'Y'))
    {
      exit(1);
    }
    fclose(fpout_ptr);
 }
*/
 
/*-----Read input data------------------------------------------------*/
 
 //fpin_ptr   = fopen(argv[1],"r");
/*
 if (fpin_ptr == NULL)
 {
    printf("\nCannot open input file\n");
    exit(1);
 }*/
  
 //fscanf(fpin_ptr,"%d",&strategy);       /*---choice of strategy-----------------*/
 //fscanf(fpin_ptr,"%d",&genmax);         /*---maximum number of generations------*/
 //fscanf(fpin_ptr,"%d",&refresh);        /*---output refresh cycle---------------*/
 //fscanf(fpin_ptr,"%d",&D);              /*---number of parameters---------------*/
 //fscanf(fpin_ptr,"%d",&NP);             /*---population size.-------------------*/
 //fscanf(fpin_ptr,"%lf",&inibound_h);    /*---upper parameter bound for init-----*/
 //fscanf(fpin_ptr,"%lf",&inibound_l);    /*---lower parameter bound for init-----*/
 //fscanf(fpin_ptr,"%lf",&F);             /*---weight factor----------------------*/
 //fscanf(fpin_ptr,"%lf",&CR);            /*---crossing over factor---------------*/
 //fscanf(fpin_ptr,"%d",&seed);           /*---random seed------------------------*/
// 目前已經採用 strategy 3 可以得到最佳結果
  strategy = 3;
  genmax = 2000;
  refresh = 100;
  // 配合機構尺寸合成, 每一個體有 9 個機構尺寸值與 5 個通過點角度值
  D = 19;
  NP = 200;
  inibound_h = 50.;
  inibound_l = 0.;
/*得到最佳解
  F = 0.85;
CR 必須介於 0 to 1. 之間
  CR = 1.;
*/
  F = 0.85;
  CR = 1.;
  seed = 3;
 
 //fclose(fpin_ptr);
 
/*-----Checking input variables for proper range----------------------------*/
 
  if (D > MAXDIM)
  {
     printf("\nError! D=%d > MAXDIM=%d\n",D,MAXDIM);
     exit(1);
  }
  if (D <= 0)
  {
     printf("\nError! D=%d, should be > 0\n",D);
     exit(1);
  }
  if (NP > MAXPOP)
  {
     printf("\nError! NP=%d > MAXPOP=%d\n",NP,MAXPOP);
     exit(1);
  }
  if (NP <= 0)
  {
     printf("\nError! NP=%d, should be > 0\n",NP);
     exit(1);
  }
  if ((CR < 0) || (CR > 1.0))
  {
     printf("\nError! CR=%f, should be ex [0,1]\n",CR);
     exit(1);
  }
  if (seed <= 0)
  {
     printf("\nError! seed=%d, should be > 0\n",seed);
     exit(1);
  }
  if (refresh <= 0)
  {
     printf("\nError! refresh=%d, should be > 0\n",refresh);
     exit(1);
  }
  if (genmax <= 0)
  {
     printf("\nError! genmax=%d, should be > 0\n",genmax);
     exit(1);
  }
  if ((strategy < 0) || (strategy > 10))
  {
     printf("\nError! strategy=%d, should be ex {1,2,3,4,5,6,7,8,9,10}\n",strategy);
     exit(1);
  }
  if (inibound_h < inibound_l)
  {
     printf("\nError! inibound_h=%f < inibound_l=%f\n",inibound_h, inibound_l);
     exit(1);
  }
 
 
/*-----Open output file-----------------------------------------------*/
 
   //fpout_ptr   = fopen(argv[2],"w");  /* open output file for writing */
 
   //if (fpout_ptr == NULL)
   //{
      //printf("\nCannot open output file\n");
      //exit(1);
   //}
 
 
/*-----Initialize random number generator-----------------------------*/
 
 rnd_uni_init = -(long)seed;  /* initialization of rnd_uni() */
 nfeval       =  0;  /* reset number of function evaluations */
 
 
 
/*------Initialization------------------------------------------------*/
/*------Right now this part is kept fairly simple and just generates--*/
/*------random numbers in the range [-initfac, +initfac]. You might---*/
/*------want to extend the init part such that you can initialize-----*/
/*------each parameter separately.------------------------------------*/
 
   for (i=0; i<NP; i++)
   {
      for (j=0; j<D; j++) /* spread initial population members */
      {
        c[i][j] = inibound_l + rnd_uni(&rnd_uni_init)*(inibound_h - inibound_l);
      }
      cost[i] = evaluate(D,c[i],&nfeval); /* obj. funct. value */
   }
   cmin = cost[0];
   imin = 0;
   for (i=1; i<NP; i++)
   {
     if(MAXIMAPROBLEM == 1)
     {
       // 改為最大化
        if (cost[i]>cmin)
        {
          cmin = cost[i];
          imin = i;
        }
      }
      else
      {
        // 最小化問題
        if (cost[i]<cmin)
        {
          cmin = cost[i];
          imin = i;
        }
      }
   }
 
   assignd(D,best,c[imin]);            /* save best member ever          */
   assignd(D,bestit,c[imin]);          /* save best member of generation */
 
   pold = &c; /* old population (generation G)   */
   pnew = &d; /* new population (generation G+1) */
 
/*=======================================================================*/
/*=========Iteration loop================================================*/
/*=======================================================================*/
 
   gen = 0;                          /* generation counter reset */
   while ((gen < genmax) /*&& (kbhit() == 0)*/) /* remove comments if conio.h */
   {                                            /* is accepted by compiler    */
      gen++;
      imin = 0;
 
      for (i=0; i<NP; i++)         /* Start of loop through ensemble  */
      {
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 2 !!!     */
       r1 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while(r1==i);            
 
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 3 !!!     */
       r2 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r2==i) || (r2==r1));
 
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 4 !!!     */
       r3 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r3==i) || (r3==r1) || (r3==r2));
 
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 5 !!!     */
       r4 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));
 
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 6 !!!     */
       r5 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));
 
 
/*=======Choice of strategy===============================================================*/
/*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
/*=======DE :  stands for Differential Evolution==========================================*/
/*=======x  :  a string which denotes the vector to be perturbed==========================*/
/*=======y  :  number of difference vectors taken for perturbation of x===================*/
/*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
/*                                                                                        */
/*=======There are some simple rules which are worth following:===========================*/
/*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
/*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
/*=======3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=*/
/*           happens.                                                                     */
/*=======4)  If you increase NP, F usually has to be decreased============================*/
/*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/
 
 
/*=======EXPONENTIAL CROSSOVER============================================================*/
 
/*-------DE/best/1/exp--------------------------------------------------------------------*/
/*-------Our oldest strategy but still not bad. However, we have found several------------*/
/*-------optimization problems where misconvergence occurs.-------------------------------*/
     if (strategy == 1) /* strategy DE0 (not in our paper) */
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D);
       L = 0;
       do
       {                       
         tmp[n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand/1/exp-------------------------------------------------------------------*/
/*-------This is one of my favourite strategies. It works especially well when the-------*/
/*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
/*-------as a first guess.---------------------------------------------------------------*/
     else if (strategy == 2) /* strategy DE1 in the techreport */
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D);
       L = 0;
       do
       {                       
         tmp[n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
/*-------should play around with all three control variables.----------------------------*/
     else if (strategy == 3) /* similiar to DE2 but generally better */
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                       
         tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
     else if (strategy == 4)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                           
         tmp[n] = bestit[n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
     else if (strategy == 5)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                           
         tmp[n] = (*pold)[r5][n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
 
/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/
 
/*-------DE/best/1/bin--------------------------------------------------------------------*/
     else if (strategy == 6) 
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand/1/bin-------------------------------------------------------------------*/
     else if (strategy == 7) 
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
     else if (strategy == 8) 
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/best/2/bin--------------------------------------------------------------------*/
     else if (strategy == 9)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = bestit[n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand/2/bin--------------------------------------------------------------------*/
     else
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = (*pold)[r5][n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         }
         n = (n+1)%D;
           }
     }
 
 
/*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/
 
     trial_cost = evaluate(D,tmp,&nfeval);  /* Evaluate new vector in tmp[] */
   if(MAXIMAPROBLEM == 1)
   {
    // 改為最大化
       if (trial_cost >= cost[i])   /* improved objective function value ? */
       {                                  
          cost[i]=trial_cost;         
          assignd(D,(*pnew)[i],tmp);
          if (trial_cost>cmin)          /* Was this a new minimum? */
          {                               /* if so...*/
             cmin=trial_cost;           /* reset cmin to new low...*/
             imin=i;
             assignd(D,best,tmp);           
          }                           
       }                            
       else
       {
          assignd(D,(*pnew)[i],(*pold)[i]); /* replace target with old value */
       }
    }
    else
    {
          // 最小化問題
       if (trial_cost <= cost[i])   /* improved objective function value ? */
       {                                  
          cost[i]=trial_cost;         
          assignd(D,(*pnew)[i],tmp);
          if (trial_cost<cmin)          /* Was this a new minimum? */
          {                               /* if so...*/
             cmin=trial_cost;           /* reset cmin to new low...*/
             imin=i;
             assignd(D,best,tmp);           
          }                           
       }                            
       else
       {
          assignd(D,(*pnew)[i],(*pold)[i]); /* replace target with old value */
       }
    }
     
      }   /* End mutation loop through pop. */
                        
      assignd(D,bestit,best);  /* Save best population member of current iteration */
 
      /* swap population arrays. New generation becomes old one */
 
      pswap = pold;
      pold  = pnew;
      pnew  = pswap;
 
/*----Compute the energy variance (just for monitoring purposes)-----------*/
 
      cmean = 0.;          /* compute the mean value first */
      for (j=0; j<NP; j++)
      {
         cmean += cost[j];
      }
      cmean = cmean/NP;
 
      cvar = 0.;           /* now the variance              */
      for (j=0; j<NP; j++)
      {
         cvar += (cost[j] - cmean)*(cost[j] - cmean);
      }
      cvar = cvar/(NP-1);
 
 
/*----Output part----------------------------------------------------------*/
 
      if (gen%refresh==1)   /* display after every refresh generations */
      { /* ABORT works only if conio.h is accepted by your compiler */
    printf("\n\n                         PRESS ANY KEY TO ABORT"); 
    printf("\n\n\n Best-so-far cost funct. value=%-15.10g\n",cmin);
 
    for (j=0;j<D;j++)
    {
      printf("\n best[%d]=%-15.10g",j,best[j]);
    }
    printf("\n\n Generation=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
    printf("\n NP=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.5g\n",
               NP,F,CR,cvar);
      }
 
      fprintf(fpout_ptr,"%ld   %-15.10g\n",nfeval,cmin);
   }
/*=======================================================================*/
/*=========End of iteration loop=========================================*/
/*=======================================================================*/
 
/*-------Final output in file-------------------------------------------*/
 
 
   fprintf(fpout_ptr,"\n\n\n Best-so-far obj. funct. value = %-15.10g\n",cmin);
 
   for (j=0;j<D;j++)
   {
     fprintf(fpout_ptr,"\n best[%d]=%-15.10g",j,best[j]);
   }
   fprintf(fpout_ptr,"\n\n Generation=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
   fprintf(fpout_ptr,"\n NP=%d    F=%-4.2g    CR=%-4.2g    cost-variance=%-10.5g\n",
           NP,F,CR,cvar); 
 
  fclose(fpout_ptr);
    
  /* Code you want timed here */
  printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
   return(0);
}
 
/*-----------End of main()------------------------------------------*/
 
// 適應函式 fittness function (cost function)
double evaluate(int D, double tmp[], long *nfeval)
{
  // 先處理通過 5 個點的四連桿問題
  // x0, y0 為左方固定點座標, 必須在 0 ~ 100 間 - 設為 tmp[0], tmp[1]
  // x1, y1 為右方固定點座標, 必須在 0 ~ 100 間 - 設為 tmp[2], tmp[3]
  // L1 為第一桿件的長度, 必須 > 0, 且小於 100 - 設為 tmp[4]
  // L2 為第二桿件的長度, 必須 > 0, 且小於 100 - 設為 tmp[5]
  // L3 為第三桿件的長度, 必須 > 0, 且小於 100 - 設為 tmp[6]
  // L5, L6 與 L2 共同圍出可動桿 L2 對應的三角形, 關注的三角形頂點即 L5 與 L6 的交點, 而 angle3 則為 L6 之對應角(為固定值)
  // L5, L6 必須 > 0, 且小於 100 - 設為 tmp[7], tmp[8]
  // 以下的角度輸入值, 會隨著目標點數的增加而增加, 其索引值由 9 + 通過點數 - 1 決定, 5 點, 則索引至 13, 若通過 25 點, 則索引值為 9 + 24 = 33
  // input_angles[] 為五的輸入的雙浮點數角度值,代表個體的角度向量值 - 分別設為 tmp[9], tmp[10], tmp[11], tmp[12], tmp[13]
  // output_points[] 則為與 input_angles[] 對應的五個三角形的頂點值, 為座標結構值, 分別有 x 與 y 分量值
  // 當利用個體的向量值, 代入 mechanism 後所得到得 output_points[] 再與 target_points[] 進行 cost function 的誤差值最小化
  /* void mechanism(double x0, double y0, double x1, double y1, double L1,
  double L2, double L3, double L5, double L6, double input_angles[NUM_OF_POINTS], struct Coord output_points[NUM_OF_POINTS])*/
  struct Coord target_points[NUM_OF_POINTS], output_points[NUM_OF_POINTS];
  double input_angles[NUM_OF_POINTS], result;
  int i;
   
  (*nfeval)++;
   
  target_points[0].x = 1.0;
  target_points[0].y = 1.0;
   
  target_points[1].x = 2.0;
  target_points[1].y = 2.0;
   
  target_points[2].x = 3.0;
  target_points[2].y = 3.0;
   
  target_points[3].x = 4.0;
  target_points[3].y = 4.0;
   
  target_points[4].x = 5.0;
  target_points[4].y = 5.0;
   
  target_points[5].x = 6.0;
  target_points[5].y = 6.0;
   
  target_points[6].x = 7.0;
  target_points[6].y = 7.0;
   
  target_points[7].x = 8.0;
  target_points[7].y = 8.0;
   
  target_points[8].x = 9.0;
  target_points[8].y = 9.0;
   
  target_points[9].x = 10.0;
  target_points[9].y = 10.0;
   
  // 輸入角度值與 tmp[] 的設定
  for(i = 0; i < NUM_OF_POINTS; i++)
  {
    input_angles[i] = tmp[i + 9];
  }
  // 呼叫 mechanism() 以便計算 output_points[]
  mechanism(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], input_angles, output_points);
   
  // for debug
  /*
  if(*nfeval%3000 == 0)
  {
    for(i = 0; i < NUM_OF_POINTS; i++)
    {
      printf("%-15.10g : %-15.10g\n", output_points[i].x, output_points[i].y);
    }
    printf("#####################################\n");
  }
  */
  // 利用 output_points[] 與 target_points 計算誤差值, 該誤差值就是 cost
  result = error_function(output_points, target_points);
  // 這裡要分別針對各變數的約束條件與範圍值來設定傳回 PENALITY 或 誤差值 result
   
  // x0 與 x1 點位於 -500 與 500 中間
    for(i = 0; i < 4; i++)
  {
    if(tmp[i] < -50 || tmp[i] > 50){
      return PENALITY;
    }
  }
   
  // 三個連桿值, 一定要為正
    for(i = 4; i < 7; i++)
  {
    if(tmp[i] < 0 || tmp[i] > 50){
      return PENALITY;
    }
  }
   
    // L5 L6 可以為 0 或負值
    for(i = 7; i < 9; i++)
  {
    if(tmp[i] < -50 || tmp[i] > 50){
      return PENALITY;
    }
  }
   
  // 角度值一定要大於 0
 
  for(i = 1; i <= NUM_OF_POINTS; i++)
  {
    if((tmp[D-i] < 0)){
      return PENALITY;
    }
  }
   
  return result;
 
  /*
   double result=0, surface = 80.0, z, volume, penality;
   (*nfeval)++;
   z = (surface-tmp[0]*tmp[1])/(2.0*(tmp[0]+tmp[1]));
   volume = tmp[0]*tmp[1]*z;
   
  if(volume <= 0){
    return PENALITY;
  }
   
  if((tmp[0] <= inibound_l)|| (tmp[0] >inibound_h)){
    return PENALITY;
  }
   
  if((tmp[1] <= inibound_l) || (tmp[1] >inibound_h)){
    return PENALITY;
  }
  // volume must >0 and max volume
  // 目前為最小化問題
   return 1+1/(volume*volume);
   */
}
 
struct Coord triangletip_coord( double x0, double y0, double R0, double R1, double x1, double y1, double localt)
{
    struct Coord tip_coord;
 
    if (localt>=0 && localt <PI)
    {
        // 目前蓋掉的式子為利用手動代換出來的版本
        //x_value = ((x1-x0)*(-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2))/(2*sqrt(pow((y1-y0),2)+pow((x1-x0),2)))-(y1-y0)*(-mech_loop)*sqrt(fabs(pow(R0,2)-pow((-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2)),2)/(4*(pow((y1-y0),2)+pow((x1-x0),2))))))/sqrt(pow((y1-y0),2)+pow((x1-x0),2))+x0;
        // 以下的式子,先利用文字編輯器,將原先 stringout() 出來的  sqrt 替換成  sqrtt, 以防止被  maxima 中的 subst("^"=pow,expr) 所替換, subst 之後,再使用文字編輯器換回來,就可以得到正確的 C 對應碼.
        tip_coord.x = pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*(mech_loop*(y1-y0)*sqrt(fabs(pow(R0,2)-pow(pow(y1-y0,2)+
    pow(x1-x0,2),-1)*pow(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2),2)/4))+(x1-x0)*pow(sqrt(pow(y1-y0,2)+
    pow(x1-x0,2)),-1)*(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2))/2)+x0;
    }
    else
    {
        // 目前蓋掉的式子為利用手動代換出來的版本
        //x_value = ((x1-x0)*(-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2))/(2*sqrt(pow((y1-y0),2)+pow((x1-x0),2)))-(y1-y0)*(mech_loop)*sqrt(fabs(pow(R0,2)-pow((-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2)),2)/(4*(pow((y1-y0),2)+pow((x1-x0),2))))))/sqrt(pow((y1-y0),2)+pow((x1-x0),2))+x0;
        tip_coord.x = pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*(-mech_loop*(y1-y0)*sqrt(fabs(pow(R0,2)-pow(pow(y1-y0,2)+
    pow(x1-x0,2),-1)*pow(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2),2)/4))+(x1-x0)*pow(sqrt(pow(y1-y0,2)+
    pow(x1-x0,2)),-1)*(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+
    pow(x1-x0,2))/2)+x0;
    }
 
// 請注意,與 Maxma 公式中的差異為,在 sqrt()中加入 fabs(),避免因為sqrt()中的負值而造成 NaN (Not a Number 問題.
    if (localt>=0 && localt <PI)
    {
        tip_coord.y = /*((x1-x0)*(-mech_loop)*sqrt(
                fabs(pow(R0,2)-pow((-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2)),2)
                /(4*(pow((y1-y0),2)+pow((x1-x0),2)))
                ))
                +(y1-y0)*(-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2))/(2*sqrt(pow((y1-y0),2)+pow((x1-x0),2))))/sqrt(pow((y1-y0),2)+pow((x1-x0),2))
                +y0;*/
                // 利用 sqrtt 居中進行代換所得到的式子
                pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*((y1-y0)*pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*(-pow(R1,2)+
    pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2))/2-mech_loop*(x1-x0)*sqrt(fabs(pow(R0,2)-pow(pow(y1-y0,2)+
    pow(x1-x0,2),-1)*pow(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2),2)/4)))+y0;
 
    }
    else
    {
        tip_coord.y = /*((x1-x0)*(mech_loop)*sqrt(
                fabs(pow(R0,2)-pow((-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2)),2)
                /(4*(pow((y1-y0),2)+pow((x1-x0),2)))
                ))
                +(y1-y0)*(-pow(R1,2)+pow(R0,2)+pow((y1-y0),2)+pow((x1-x0),2))/(2*sqrt(pow((y1-y0),2)+pow((x1-x0),2))))/sqrt(pow((y1-y0),2)+pow((x1-x0),2))
                +y0;*/
                pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*((y1-y0)*pow(sqrt(pow(y1-y0,2)+pow(x1-x0,2)),-1)*(-pow(R1,2)+
    pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2))/2+mech_loop*(x1-x0)*sqrt(fabs(pow(R0,2)-pow(pow(y1-y0,2)+
    pow(x1-x0,2),-1)*pow(-pow(R1,2)+pow(R0,2)+pow(y1-y0,2)+pow(x1-x0,2),2)/4)))+y0;
    }
   
  return tip_coord;
}
 
double distance(double x0, double y0, double x1, double y1)
{
    double distance_value;
    distance_value = sqrt(pow((x1-x0),2) + pow((y1-y0),2));
    return distance_value;
}
 
double rr(double L1, double dd, double theta)
{
    double rr_value;
    rr_value = sqrt(L1*L1+dd*dd-2*L1*dd*cos(theta));
    return rr_value;
}
 
// 輸入每一個體的變數向量, 然後求各三角形頂點的座標陣列[NUM_OF_POINTS]
void mechanism(double x0, double y0, double x1, double y1, double L1,
  double L2, double L3, double L5, double L6, double input_angles[NUM_OF_POINTS], struct Coord output_points[NUM_OF_POINTS])
{
  // 此函式要輸入控制變數, 然後計算機構尺寸合成的關注點座標位置
  // 以下為可能的處理變數宣告
  // 這裡希望能夠定義一個 struct 來處理座標點
  double rr_length, dd_length, angle;
  struct Coord link1_tip, link2_tip, triangle_tip;
    double angle2, angle3;
  int i;
   
  // 開始進行三角形頂點座標的計算
  // 以下變數由每一個體向量提供
  /*
    x0 = 0.0;
    y0 = 0.0;
    x1 = 10.0;
    y1 = 0.0;
    L1 = 5.0;
    L2 = 20;
    L3 = 10;
    L5 = 10;
    L6 = 10;
  */
  dd_length = distance(x0, y0, x1, y1);
  /* 設法表示 triangle 所對應的 local 角度,表示為已知變數與 t 的函式 */
  angle3 = acos((pow(L2,2)+pow(L5,2)-pow(L6,2))/(2*L2*L5));
   
  for(i = 0; i < NUM_OF_POINTS; i++)
  {
    // 先建立第一點座標, 即 i=0 者
    // i=0;
    // angle = i*degree;
    /*
    // 利用角度增量進行運算, 相對於 input_angles[0] 作為基準
    if(i > 0)
    {
      input_angles[i] = input_angles[i] + input_angles[i-1];
    }
    */
    angle = input_angles[i]*degree;
    rr_length = rr(L1, dd_length, angle);
    // 第一次三角形疊代
    link1_tip = triangletip_coord(x0, y0, L1, rr_length, x1, y1, angle);
    // 第二次三角形疊代
    /* 設法表示 link2 所對應的 local 角度,表示為已知變數與 t 的函式 */
    angle2 = acos((pow(L2,2)+pow(rr_length,2)-pow(L3,2))/(2*L2*rr_length));
    link2_tip = triangletip_coord(link1_tip.x, link1_tip.y, L2, L3, x1, y1, angle2);
    // 第三次三角形疊代
    //triangle_tip = triangletip_coord(link1_tip.x, link1_tip.y, L5, L6, link2_tip.x, link2_tip.y, angle3);
    // output_points[i] = triangletip_coord(link1_tip.x, link1_tip.y, L5, L6, link2_tip.x, link2_tip.y, angle3);
    // 這裡要嘗試利用 finaltip_coord() 求 tip3 座標, 而 L5 與 L6 可 0 可負
    output_points[i] = finaltip_coord(link1_tip, link2_tip, L5, L6);
  }
}
 
double error_function(struct Coord output_points[NUM_OF_POINTS], struct Coord target_points[NUM_OF_POINTS])
{
  double error = 0.0;
  int i;
  for(i = 0; i < NUM_OF_POINTS; i++)
  {
    error += fabs(distance(output_points[i].x, output_points[i].y, target_points[i].x, target_points[i].y));
  }
  return error;
}
 
struct Coord finaltip_coord(struct Coord tip1_coord, struct Coord tip2_coord, double r1, double r2)
{
  struct Coord tip3_coord;
  double theta3, theta4, length3, length4;
  length3 = sqrt(pow(tip2_coord.x - tip1_coord.x,2) + pow(tip2_coord.y - tip1_coord.y,2));
  length4 = sqrt(pow(r1,2) + pow(r2,2));  
  theta3 = acos((tip2_coord.x - tip1_coord.x) / length3);
  theta4 = acos(r1 / length4);
  tip3_coord.x = tip1_coord.x + length4 * cos(theta3 + theta4);
  tip3_coord.y = tip1_coord.y + length4 * sin(theta3 + theta4);
   
  return tip3_coord;
}
