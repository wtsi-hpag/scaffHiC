/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2018  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of ScaffHiC pipeline.                                 *
 *                                                                          *
 *  ScaffHiC is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

#define PADCHAR '-'
#define Max_N_NameBase 60
static char **S_Name,**R_Name;
static int *hit_rddex,*hit_score,*hit_rcdex,*hit_locus1,*hit_locus2,*hit_matlocu1,*hit_matlocu2,*hit_matindex;
static int *ctg_length,*hit_index;

/* SSAS default parameters   */
static int IMOD=0;
static int nContig=0;
static int n_lattice = 10000;
static int belt_size = 2;
static int n_file = 0;
static int max_len = 0;
static float m_score = 400.0;
static int n_contigs = 2;
static long G_Size = 0;

int main(int argc, char **argv)
{
    FILE *namef;
    int i,nSeq,args,idt;
    int n_contig,n_reads,nseq;
    void Matrix_Process(char **argv,int args,int nSeq);
    char *st,*ed;
    char line[2000]={0},tempc1[60],rdname[60];
    char **cmatrix(long nrl,long nrh,long ncl,long nch);

    if(argc < 2)
    {
      printf("Usage: %s -belt 2 -lattice 10000  <Input_Mapping_file> <Output_file1> <Output_file2>\n",argv[0]);

      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&n_file);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-belt"))
       {
         sscanf(argv[++i],"%d",&belt_size);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-lattice"))
       {
         sscanf(argv[++i],"%d",&n_lattice);
         args=args+2;
       }
    }

    nseq=0;
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 

    n_contigs = nseq;
    if((ctg_length = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - scf_list\n");
      exit(1);
    }

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
/*  read the alignment files         */
    i=0;
    G_Size = 0;
    while(fscanf(namef,"%s %s %d %s",tempc1,tempc1,&ctg_length[i],tempc1)!=EOF)
    {
        if(ctg_length[i] > max_len)
          max_len = ctg_length[i];
        G_Size = G_Size+ctg_length[i];
        i++;
    }
    fclose(namef);

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 

    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus1\n");
      exit(1);
    }
    if((hit_locus2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((hit_index = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_index\n");
      exit(1);
    }

    nSeq=nseq;
    n_contig=0;
    n_reads=0;

    printf("reads: %d %s\n",nseq,argv[args]);
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%d %d %d",&hit_index[i],&hit_locus1[i],&hit_locus2[i])!=EOF)
    {
        i++;
    }
    fclose(namef);

    n_reads=i;


    printf("read2: %d %s\n",i,argv[args]);
    Matrix_Process(argv,args,n_reads);

    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Matrix_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,nseq,n_blocks,stopflag;
     int idt,idt1,n_zeros,n_cells;
     FILE *namef,*namef2;
     char *st,*ed,tempc1[60],rdname[60],rdname1[60],ctgname[60],line[2000];
     int **heat_map,**s_matrix,*ctg_locus,*ctg_index,*ctg_rddex;
     int *ctg_lnzero,*ctg_btzero,*ctg_dbzero;
     int num_hits,num_hit1,num_hit2,rcdex,rsize,rsize2,size_row;
     int **imatrix(long nrl,long nrh,long ncl,long nch);
     int *cover_genome,*ctg_idex1,*ctg_idex2,*ctg_mask,*ctg_rcdex1;
     void ArraySort2_Int2(int n, int *arr, int *brr);
     double sqrt(double arg);
     float rate;

     size_row = 20;

/*
     nseq=0;
     if((namef = fopen(argv[args+2],"r")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     while(!feof(namef))
     {
       if(fgets(line,2000,namef) == NULL)
       {
       }
       if(feof(namef)) break;
       nseq++;
     }
     fclose(namef);

     if((ctg_locus = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_locus\n");
       exit(1);
     }
     if((ctg_index = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_index\n");
       exit(1);
     }
     if((ctg_rddex = (int *)calloc(n_contigs,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_rddex\n");
       exit(1);
     }
     s_matrix=imatrix(0,n_contigs,0,n_contigs);

     if((namef = fopen(argv[args+2],"r")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }

     i=0;
     idt = 0;
     idt1 = 0;
     while(fscanf(namef,"%s %s %d %s %s %s",rdname,ctgname,&ctg_locus[i],tempc1,tempc1,tempc1)!=EOF)
     {
         st = ctgname;
         ed = strrchr(ctgname,'_');
         idt = atoi(ed+1);
         ctg_index[i] = atoi(ed+1);
         if(strcmp(rdname1,rdname) == 0)
         {
           s_matrix[idt][idt1]++;
           s_matrix[idt1][idt]++;
         }
         strcpy(rdname1,rdname);
         idt1 = idt;
         i++;
     }
     fclose(namef);

     for(i=0;i<n_contigs;i++)
     {
        for(j=0;j<n_contigs;j++)
           ctg_rddex[j] = j;
        ArraySort2_Int2(n_contigs,s_matrix[i],ctg_rddex);
        printf("Links: %d %d",i,ctg_length[i]);
        for(k=0;k<size_row;k++)
           printf(" %d ",s_matrix[i][k]);
        printf("\n");
     }

                     */
     rsize = 100000; 
     heat_map = imatrix(0,rsize,0,rsize);
     n_blocks = rsize*rsize; 

     n_cells = (max_len+40000)/20000; 
          printf("cell: %d\n",n_cells);
     if((ctg_lnzero = (int *)calloc(n_cells,sizeof(int))) == NULL)
     {
        printf("fmate: calloc - ctg_index\n");
        exit(1);
     }
     if((ctg_btzero = (int *)calloc(n_cells,sizeof(int))) == NULL)
     {
        printf("fmate: calloc - ctg_index\n");
        exit(1);
     }
     if((ctg_dbzero = (int *)calloc(n_cells,sizeof(int))) == NULL)
     {
        printf("fmate: calloc - ctg_index\n");
        exit(1);
     }
     for(i=0;i<nSeq;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(hit_index[i]==hit_index[j])
          {
            j++;
          }
          else
            stopflag=1;
        }
        num_hits = j-i;
          printf("chr: %d\n",num_hits);
        if(num_hits>=10) 
        {
          int chr_len = ctg_length[hit_index[i]];

          rate = num_hits;
          rate = rate/chr_len;

          if(rate <= 0.02)
            n_lattice = 30000;
          else if((rate > 0.02)&&(rate < 0.05))
            n_lattice = 25000;
          else if((rate >= 0.05)&&(rate < 0.1))
          {
            float rate2 = 0.05;
            rate2 = sqrt(rate2/rate)*(rate2/rate)*15000;
            n_lattice = rate2/1000;
            n_lattice = n_lattice*1000; 
          }
          else if((rate >= 0.1)&&(rate <0.2))
            n_lattice = 4000;
          else if(rate >= 0.2)
            n_lattice = 3000;

          n_lattice = 20000;
          n_cells = ctg_length[hit_index[i]]/n_lattice;


          for(k=0;k<n_cells;k++)
          {
             ctg_lnzero[k] = 0;
             ctg_btzero[k] = 0;
             for(n=0;n<n_cells;n++)
                heat_map[k][n] = 0;;
          }
          for(n=i;n<j;n++)
          {
             int idi,idj;
             idi = hit_locus1[n]/n_lattice;
             idj = hit_locus2[n]/n_lattice;
             heat_map[idi][idj]++;
             heat_map[idj][idi]++;
          }
          printf("Scaffold: %d %d %d %f %d\n",hit_index[i],belt_size,chr_len,rate,n_lattice);
          for(k=0;k<n_cells;k++)
          {
             int sum_gap = 0;
            
             if(k < belt_size)
             {
               for(n=0;n<(k-2);n++)
                  sum_gap = sum_gap + heat_map[k][n];
             }
             else
             {
               for(n=(k-belt_size);n<(k-10);n++)
                  sum_gap = sum_gap + heat_map[k][n];
             }
             if(n_file == 0)
             {
                int end_cell = 0;
                int end_loci = 0;
                int end_cell2 = 0;
                int end_loci2 = 0;
                int num_sums;
                int line_loci,line_cell,line_loci2,line_cell2,row_loci,row_cell;
                int row_loci2,row_cell2,z_loci,z_cell,z_loci2,z_cell2;
                int l_zeros,r_zeros,l10_zeros,r10_zeros,l20_zeros,d20_zeros;

                for(n=0;n<n_cells;n++)
                {
                   if(heat_map[k][n] >= 1)
                   {
                     end_cell = heat_map[k][n];
                     end_loci = n;
                     if(heat_map[k][n] >= 2)
                     {
                       end_cell2 = heat_map[k][n];
                       end_loci2 = n;
                     } 
                   } 
                }
                line_loci = end_loci;
                line_cell = end_cell;
                line_loci2 = end_loci2;
                line_cell2 = end_cell2;

                num_sums = 0;
                for(n=(k+3);n<n_cells;n++)
                   num_sums = num_sums+heat_map[k][n];
                
                num_sums = num_sums/n_cells;

                end_cell = 0;
                end_loci = n_cells;
                n_zeros = 0;
                for(n=0;n<(n_cells-5);n++)
                {
                   if((heat_map[k][n] == 0)&&(heat_map[k][n+1] == 0)&&(heat_map[k][n+2] == 0)&&(heat_map[k][n+3] == 0)&&(heat_map[k][n+4] == 0))
                   {
                     end_cell = 5;
                     end_loci = n;
                     n_zeros++;
                   } 
                }
                z_loci = end_loci;
                z_cell = end_cell;

                l_zeros = 0;
                l10_zeros = 0;
                l20_zeros = 0;
                d20_zeros = 0;


                for(m=(k-5);m<(k+5);m++)
                {
                   if((m>=5)&&(m<(n_cells-6))&&(heat_map[m+5][m-5] == 0))
                     d20_zeros++;
                }
                for(n=(k+3);n<(n_cells-1);n++)
                {
                   if((n<(k+13))&&(heat_map[n][k] == 0))
                     l10_zeros++;
                   if((n<(k+23))&&(n>=(k+13))&&(heat_map[n][k] == 0))
                     l20_zeros++;
                   if(heat_map[k][n] == 0)
                     l_zeros++;
                }
 
                r_zeros = 0; 
                r10_zeros = 0; 
                for(n=(k+3);n<(n_cells-1);n++)
                {
                   if((n<(k+13))&&(heat_map[n][k] == 0))
                       r10_zeros++;
                   if(heat_map[n][k] == 0)
                     r_zeros++;
                }
                ctg_lnzero[k] = l_zeros; 
                ctg_btzero[k] = l10_zeros; 
                ctg_dbzero[k] = d20_zeros; 
                printf("Map: %d %d %d %d || %d %d %d %d || %d || %d %d %d || %d %d || %d %d %d|\n",hit_index[i],k,k*n_lattice,n_cells,line_loci,line_cell,line_loci2,line_cell2,num_sums,z_loci,z_cell,n_zeros,l_zeros,r_zeros,l10_zeros,r10_zeros,d20_zeros);
             }
             else
             {
                int sum_belt1 = 0;
                int sum_belt2= 0;
                for(n=0;n<(k-5)&&(k>=5);n++)
                   sum_belt1 = sum_belt1 + heat_map[k][n];
                for(n=(k+5);n<n_cells;n++)
                   sum_belt2 = sum_belt2 + heat_map[k][n];
                printf("Map: %d %d %d %d |%d %d| ",hit_index[i],k,k*n_lattice,sum_gap,sum_belt1,sum_belt2);
                printf("\n");
             }
/*
             for(n=0;n<n_cells;n++)
             {
                if((n!=k)&&(n<k)&&(n>(k-belt_size)))
                {
                  if((heat_map[k][n] < 10)&&((chr_len - (n*n_lattice)) > 500000)&&((n*n_lattice) > 500000))
                     printf("Map: %d %d %d %d %d %d %d\n",hit_index[i],k,n,k*n_lattice,n*n_lattice,chr_len,heat_map[k][n]);
                }
             }     */
          }
          for(k=1;k<(n_cells-1);k++)
          {
             if((ctg_btzero[k] >= 6)&&(ctg_dbzero[k] >= 4))
             {
               if((abs(ctg_lnzero[k] - ctg_lnzero[k-1]) >= 20)||(abs(ctg_lnzero[k] - ctg_lnzero[k+1]) >= 20)) 
                 printf("Break: %d %f %d %d || %d %d %d || %d %d %d || %d %d %d\n",hit_index[i],rate,k,k*n_lattice,ctg_btzero[k],ctg_btzero[k],ctg_btzero[k],ctg_lnzero[k-1],ctg_lnzero[k],ctg_lnzero[k+1],ctg_dbzero[k-1],ctg_dbzero[k],ctg_dbzero[k+1]);
             }
          }
        }
        i=j-1; 
     }
}


#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_float(int n, float *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,b,jstack=0,NSTACK=50,istack[NSTACK],MIN=7;
     float a,temp;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort_float2(int n, float *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,b,jstack=0,NSTACK=50,istack[NSTACK],MIN=7;
     float a,temp;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

/* creat char matrix with subscript ange fm[nrl...nrh][ncl...nch]  */
float   **fmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **fm;

        /* allocate pointers to rows        */
        if((fm=(float **)calloc(nrow,sizeof(float*)))==NULL)
        {
           printf("error fmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        fm+=0;
        fm-=nrl;

        /* allocate rows and set pointers to them        */
        if((fm[nrl]=(float *)calloc(nrow*ncol,sizeof(float)))==NULL)
        {
           printf("error fmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        fm[nrl]+=0;
        fm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           fm[i]=fm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return fm;
}


