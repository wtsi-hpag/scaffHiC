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
#define Max_N_NameBase2 14 
static char **R_Name;
static int *hit_rddex,*hit_score,*hit_rcdex,*hit_locus1,*superlength,*hit_matlocu1,*hit_matlocu2,*hit_matindex;
static int *ctg_length,*hit_index,*n_links,*n_list,*n_head,*map_loci;;

/* SSAS default parameters   */
static int IMOD=0;
static int nContig=0;
static int file_flag = 1;
static int ctg_minlen = 100000;
static int max_ctg = 1000000;
static int max_len = 0;
static int max_blocks = 0;
static int mtg_length = 100000000;
static int n_depth = 3;
static int c_pairs = 100;
static int breakid = 43;
static int nstep = 10000;
static int n_edges = 0;
static float m_score = 200.0;
static float d_score1 = 0.5;
static float d_score2 = 0.6;
static float c_score1 = 0.5;
static float c_score2 = 0.6;
static float all_rate = 0.0;
static int i_getindex = 3;
static int j_getindex = 52;

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
      printf("Usage: %s -depth 3 <Input_HiC-alignment> <Output_Break_Loci>\n",argv[0]);

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
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&n_depth);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&ctg_minlen);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-break"))
       {
         sscanf(argv[++i],"%d",&breakid);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%f",&m_score);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-step"))
       {
         sscanf(argv[++i],"%d",&nstep);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-index"))
       {
         sscanf(argv[++i],"%d",&i_getindex);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
    }


    fflush(stdout);
    if(system("ps aux | grep scaffHiC_break; date") == -1)
    {
//        printf("System command error:\n);
    }

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

    if((hit_rddex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - insert\n");
      exit(1);
    }
    if((hit_score = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_score\n");
      exit(1);
    }
    if((hit_rcdex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_rcdex\n");
      exit(1);
    }
    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus1\n");
      exit(1);
    }
    if((ctg_length = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_length\n");
      exit(1);
    }
    if((superlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - superlength\n");
      exit(1);
    }
    if((hit_index = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_index\n");
      exit(1);
    }
    if((map_loci = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - map_loci\n");
      exit(1);
    }

    nSeq=nseq;
    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
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
    max_ctg = 0;
    max_len = 0;
    while(fscanf(namef,"%s %s %d %d %s %d",R_Name[i],rdname,&hit_locus1[i],&hit_score[i],tempc1,&superlength[i])!=EOF)
    {
        st = rdname;
        ed = strrchr(rdname,'_');
        idt = atoi(ed+1);
        if(idt > max_ctg)
          max_ctg = idt;
        ctg_length[idt] = superlength[i];
        if(superlength[i] > max_len)
          max_len = superlength[i];
        hit_index[i] = idt;
        i++;
    }
    fclose(namef);

    max_ctg = max_ctg + 1;
    max_len = max_len + 1;
    max_blocks = max_len/nstep;
    if((n_links= (int *)calloc(max_blocks,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: n_links\n");
      exit(1);
    }
    if((n_list= (int *)calloc(max_blocks,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: n_list\n");
      exit(1);
    }
    if((n_head= (int *)calloc(max_blocks,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: n_head\n");
      exit(1);
    }

    n_reads=i;
    Matrix_Process(argv,args,n_reads);

    printf("Job finished for %d reads!\n",n_reads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Matrix_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,n_scaff,n_blocks,len_thresd;
     FILE *namef;
     long num_cells,n_Bases;
     int num_hits,num_hit1,num_hit2,rcdex,rsize,rsize2,size_row,dlinks[5];
     int stopflag,offset,*ray,*dex;
     void ArraySort_Mix(int n, long *arr, int *brr);
     void ArraySort_float2(int n, float *arr, int *brr);
     char *st,*ed;
     int **r_matrix,**rc0_matrix,**rc1_matrix,**rc2_matrix,**rc3_matrix;
     float rate;
     int **imatrix(long nrl,long nrh,long ncl,long nch);
     float **fmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     void ArraySort_Int2(int n, int *arr, int *brr);
     void Reads_Lattices(int n, int m, int *arr, int *brr);
     int *hit_linksdex,*link_locus,*link_locu2,*link_locu3,*head_locus;
     int n_length = ctg_minlen;

     rsize = max_ctg+10; 
     n_blocks = rsize*rsize;
     len_thresd = 25000000;

     if((link_locus = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_links\n");
       exit(1);
     }
     if((head_locus = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - head_locus\n");
       exit(1);
     }
     if((link_locu2 = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_links\n");
       exit(1);
     }
     if((link_locu3 = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_links\n");
       exit(1);
     }

 
     printf("contigs: %d %d %d\n",rsize,max_ctg,nSeq);
     size_row = n_depth;
     r_matrix=imatrix(0,rsize,0,rsize);
     rc0_matrix=imatrix(0,rsize,0,rsize);
     rc1_matrix=imatrix(0,rsize,0,rsize);
     rc2_matrix=imatrix(0,rsize,0,rsize);
     rc3_matrix=imatrix(0,rsize,0,rsize);
     num_hits =0;
     k = 0;
     offset = 0;
     for(i=0;i<nSeq;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(R_Name[i],R_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2) 
        {
          int idi,idt,len1,len2,loci1,loci2;

          if(hit_index[i] < hit_index[i+1])
          {
            idi = hit_index[i];
            idt = hit_index[i+1];
            len1 = superlength[i]/2;
            len2 = superlength[i+1]/2;
            loci1 = hit_locus1[i];
            loci2 = hit_locus1[i+1];
          }
          else
          {
            idi = hit_index[i+1];
            idt = hit_index[i];
            len2 = superlength[i]/2;
            len1 = superlength[i+1]/2;
            loci2 = hit_locus1[i];
            loci1 = hit_locus1[i+1];
          }

          if(loci1 < len1)
          {
            if(loci2 < len2)
            {
              rc0_matrix[idi][idt]++;
            }
            else
            {
              rc2_matrix[idi][idt]++;
            }
          }
          else
          {
            if(loci2 < len2)
            {
              rc1_matrix[idi][idt]++;
            }
            else
            {
              rc3_matrix[idi][idt]++;
            }
          } 
        }
        else
        {
          printf("www: %s %d\n",R_Name[i],superlength[i]);
        }
        i=j-1;
     }

     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
        {
           if(j>i)
           {
             rc0_matrix[j][i] = rc0_matrix[i][j];
             rc1_matrix[j][i] = rc2_matrix[i][j];
             rc2_matrix[j][i] = rc1_matrix[i][j];
             rc3_matrix[j][i] = rc3_matrix[i][j];
           }
        }
     }
     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
        {
          r_matrix[i][j] = rc0_matrix[i][j]+rc1_matrix[i][j]+rc2_matrix[i][j]+rc3_matrix[i][j];
        }
     }

     num_cells = 0;
     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
        {
           int idh = i*max_ctg;
           hit_rddex[j] = j;
           num_cells = num_cells+r_matrix[i][j];
           link_locu3[idh+j] = r_matrix[i][j];
        }
     }

     head_locus[0] = 0;
     for(i=1;i<n_blocks;i++)
     {
        head_locus[i] = head_locus[i-1] + link_locu3[i-1];
     }

     num_cells = num_cells + 20000;
     if((hit_matlocu1 = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_matlocus\n");
       exit(1);
     }
     if((hit_matlocu2 = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_matlocus\n");
       exit(1);
     }
     if((hit_matindex = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_matindex\n");
       exit(1);
     }
     if((hit_linksdex = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_linksdex\n");
       exit(1);
     }

     for(i=0;i<num_cells;i++)
        hit_linksdex[i] = i;

     for(i=0;i<nSeq;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(R_Name[i],R_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2) 
        {
          int idi,idt,ikk,len1,len2,loci1,loci2;
          int idh1,idh2;

          if(hit_index[i] < hit_index[i+1])
          {
            idi = hit_index[i];
            idt = hit_index[i+1];
            len1 = superlength[i]/2;
            len2 = superlength[i+1]/2;
            loci1 = hit_locus1[i];
            loci2 = hit_locus1[i+1];
          }
          else
          {
            idi = hit_index[i+1];
            idt = hit_index[i];
            len2 = superlength[i]/2;
            len1 = superlength[i+1]/2;
            loci2 = hit_locus1[i];
            loci1 = hit_locus1[i+1];
          }

          idh1 = idi*max_ctg+idt;
          idh2 = idt*max_ctg+idi;

          hit_matlocu1[head_locus[idh1]+link_locus[idh1]] = loci1;
          hit_matlocu1[head_locus[idh2]+link_locu2[idh2]] = loci2;
          link_locus[idh1]++;
          link_locu2[idh2]++;
        }
        else
        {
          printf("www: %s %d\n",R_Name[i],superlength[i]);
        }
        i=j-1;
     }

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
     n_Bases = 0;
     for(i=0;i<max_ctg;i++)
        n_Bases = n_Bases + ctg_length[i];

     all_rate = nSeq;
     all_rate = all_rate/n_Bases;
     printf("Total number of contigs: %d %f\n",max_ctg,all_rate);
     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
           hit_rddex[j] = j;
        ArraySort2_Int2(max_ctg,r_matrix[i],hit_rddex);
 

        printf("scaffold: %d %d %d\n",i,max_ctg,ctg_length[i]);

        printf("matrix: %d %d\n",max_ctg,ctg_length[i]);
        for(j=0;j<size_row;j++)
        {
           int idi,idj,offset1,offset2,c1_pairs,c2_pairs;
           int idd = hit_rddex[j];
           int n_pairs = r_matrix[i][j];
           int n_pairs_half = n_pairs/2;
           float rr1,rr2,rat1,rat2,sq1,sq2,rat_size;
           int halflen_1,halflen_2,half_len1,half_len2;
           int rcindex1,rcindex2,rcindex;
   
           rcindex1 = 1;
           rcindex2 = 1;

           idi = i;
           idj = hit_rddex[j]; 
           
           offset1 = head_locus[idi*max_ctg+idj];
           offset2 = head_locus[idj*max_ctg+idi];
           ray = hit_matlocu1;
           dex = hit_linksdex; 
           ArraySort_Int2(n_pairs,ray+offset1,dex+offset1);
           ArraySort_Int2(n_pairs,ray+offset2,dex+offset2);

           halflen_1 = hit_matlocu1[offset1+n_pairs_half];
           halflen_2 = hit_matlocu1[offset2+n_pairs_half];
           half_len1 = ctg_length[idi]/2;
           half_len2 = ctg_length[idj]/2;


           if(halflen_1 > ctg_length[idi]/2)
           {
             rcindex1 = 0;
             halflen_1 = ctg_length[idi]-halflen_1;
           }
           if(halflen_2 > ctg_length[idj]/2)
           {
             rcindex2 = 0;
             halflen_2 = ctg_length[idj]-halflen_2;
           }
         
           if(rcindex1 == 0)
           {
             if(rcindex2 == 0)
               rcindex = 0;
             else
               rcindex = 1;
           }
           else
           {
             if(rcindex2 == 0)
               rcindex = 2;
             else
               rcindex = 3;
           }


           rat1 = half_len1;
           rat1 = rat1/halflen_1;
           rat2 = half_len2;
           rat2 = rat2/halflen_2;

           if(rat1 > 10.0)
             rat1 = 1.0;
           if(rat2 > 10.0)
             rat2 = 1.00;

           c1_pairs = n_pairs;
           c2_pairs = n_pairs;
//           if((rat1 > 1.8)||(rat2 > 1.8)||(ctg_length[idi] < 1000000))
//           if(idi = breakid)
           i_getindex = idi;
           j_getindex = idj;
           n_edges = 0;
           d_score1 = 0;
           d_score2 = 0;  
           Reads_Lattices(n_pairs,ctg_length[idi],ray+offset1,dex+offset1);
//           if((n_edges >= 5)&&(d_score1 > c_score1)&&(d_score2 > c_score2))
//             Break_Location(n_pairs,ctg_length[idi],ray+offset1,dex+offset1);
        }
     }
     fclose(namef);
}

/* ================================================================= */
void Reads_Lattices(int nSeq, int R_len, int *rd_locus, int *rd_index)
/* ================================================================= */
{
    int i,j,k,stopflag,num_steps,num_ave,set_ave,num_ave2,set_ave2,set5;
    int num_breaks,num_hits,n_cells,max_locu,max_step,BAR = 0;
    double rate,rate2,max_rate,max_rate2,max_rate3;
    int hit_1,hit_2,hit_3,hit_4,gap_set;
    void Break_Location(int n, int m, int k, int p,int *arr, int *brr);
    FILE *namef2;

    if(R_len > 0)
    {
      rate = nSeq;
      rate = rate/R_len;
    }
    else
      rate = 0.0;
    n_cells = R_len/nstep;

    for(i=0;i<n_cells;i++)
    {
       n_links[i] = 0;
       n_list[i] = 0;
       n_head[i] = 0;
    }

    num_hits = 0;
    num_steps = 0; 
    for(i=0;i<nSeq;i++)
    {
/*     search reads with an index < i     */
/*     search reads with an index > i     */
       stopflag=0;
       j=i+1;
       while((j<nSeq)&&(stopflag==0))
       {
         if((rd_locus[j]<(BAR+nstep))&&(rd_locus[i]>=BAR)&&(rd_locus[j]<=R_len))
         {
           j++;
         }
         else
           stopflag=1;
       }
       {
         n_links[num_hits] = j-i-1;
         n_list[num_hits] = j-i-1;
         BAR = BAR+nstep;
         num_hits++;
         num_steps++;
         if(BAR >= R_len)
         {
           j = nSeq;
           i = nSeq; 
         }
       }
       i=j-1;
     }

     if(n_cells >= 1)
       num_ave = nSeq/n_cells;
     else
       num_ave = 0;

     gap_set = 0.5*num_ave;
     set_ave = 0.15*num_ave;
     num_ave2 = 0;
     num_hits = 0;
     for(i=0;i<n_cells;i++)
     {
        if(n_links[i] >= set_ave)
        {
          num_ave2 = num_ave2+n_links[i];
          num_hits++;
        } 
     }

     if(num_hits >=1)
     {
       num_ave2 = num_ave2/num_hits;
//       set_ave2 = 0.2*num_ave2;
       set_ave2 = 0.1*num_ave2;
     }
     else
     {
       num_ave2 = 0;
       set_ave2 = 0;
     }

     for(i=0;i<n_cells;i++)
     {
        if(n_links[i] < 5)
          n_links[i] = 0;
        else if(n_links[i] > (num_ave/3))
          n_links[i] = num_ave;
//        else if(n_links[i] > num_ave2)
//          n_links[i] = num_ave2;
     }
 
     for(i=0;i<n_cells;i++)
     {
        if((i>=20)&&(i<(n_cells-20))&&(n_links[i] == 0))
        {
          int ave10_1 = 0;
          int ave10_2 = 0;
          for(j=(i-16);j<=(i-6);j++)
          {
             if((j>=0)&&(n_links[j] >= gap_set))
             {
               ave10_1 = ave10_1+1;
             }
          }
          for(j=(i+6);j<(i+16);j++)
          {
             if((j < n_cells)&&(n_links[j] >= gap_set))
             {
               ave10_2 = ave10_2+1;
             }
          }
          if((ave10_1 >= 2)&&(ave10_2 >=2))
            n_links[i] = num_ave;
        }
        else if((i>=20)&&(i<(n_cells-20))&&(n_links[i] <= 10))
        {
          int ave10_1 = 0;
          int ave10_2 = 0;
          for(j=(i-16);j<=(i-6);j++)
          {
             if((j>=0)&&(n_links[j] >= gap_set))
             {
               ave10_1 = ave10_1+1;
             }
          }
          for(j=(i+6);j<(i+16);j++)
          {
             if((j < n_cells)&&(n_links[j] >= gap_set))
             {
               ave10_2 = ave10_2+1;
             }
          }
          if((ave10_1 >= 2)&&(ave10_2 >=2))
            n_links[i] = num_ave;
        }
     }
 
     for(i=0;i<n_cells;i++)
     {
        if((i>=20)&&(i<(n_cells-20))&&(n_links[i] == 0))
        {
          int ave10_1 = 0;
          int ave10_2 = 0;
          for(j=(i-16);j<=(i-6);j++)
          {
             if((j>=0)&&(n_links[j] >= gap_set))
             {
               ave10_1 = ave10_1+1;
             }
          }
          for(j=(i+6);j<(i+16);j++)
          {
             if((j < n_cells)&&(n_links[j] >= gap_set))
             {
               ave10_2 = ave10_2+1;
             }
          }
          if((ave10_1 >= 2)&&(ave10_2 >=2))
            n_links[i] = num_ave;
        }
     }
 
     set5 = 0;
     for(i=0;i<10;i++)
        set5 = set5+n_links[i];
     set5 = set5/10;


     printf("Set: %d %d %d || %f || %d %d\n",set5,set_ave,set_ave2,rate,num_ave,num_ave2);

     max_rate = 0.0;
     max_rate2 = 0.0;
     max_rate3 = 0.0;
     max_locu = 0;
     max_step = 0;

     hit_1 = 0;
     hit_2 = 0;
     hit_3 = 0;
     hit_4 = 0;
     for(i=0;i<n_cells;i++)
     {
        int ave5_1 = 0;
        int ave5_2 = 0;
        int cutave_1 = 0;
        int cutave_2 = 0;
        float rr = 0.0;
        float rm = 0.0;
        float rn = 0.0;

        for(j=0;j<i;j++)
        {
           if(j>=0)
             cutave_1 = cutave_1 + n_links[j];
        }
        if(i>=1)
          cutave_1 = cutave_1/i;
        else
          cutave_1 = 0;

        for(j=(i+1);j<n_cells;j++)
        {
           cutave_2 = cutave_2 + n_links[j];
        }
        if((n_cells-i)>=1)
          cutave_2 = cutave_2/(n_cells-i);
        else
          cutave_2 = 0;

//        for(j=(i-6);j<=(i-2);j++)
        for(j=(i-11);j<=(i-2);j++)
        {
           if(j>=0)
             ave5_1 = ave5_1 + n_links[j];
        }
        ave5_1 = ave5_1/10; 
        for(j=(i+2);j<(i+12);j++)
        {
           if(j<n_cells)
             ave5_2 = ave5_2 + n_links[j];
        }
        ave5_2 = ave5_2/10;
        if(ave5_1 > ave5_2)
        {
          rr =  ave5_1 - ave5_2;
          rm = rr/ave5_1;
//          rn = rr/(ave5_1+num_ave);
          rn = rr/(rr+num_ave*0.5);
        }
        else
        {
          rr =  ave5_2 - ave5_1;
          rm = rr/ave5_2;
//          rn = rr/(ave5_2+num_ave);
          rn = rr/(rr+num_ave*0.5);
        }
        rr = rr/num_ave;
        if((rn > max_rate)&&(i>=12)&&(i<(n_cells-12)))
        {
          max_rate = rn;
          max_locu = i*nstep;
          max_step = i;
          max_rate2 = rm;
          max_rate3 = rr;
          hit_1 = ave5_1;
          hit_2 = ave5_2;
          hit_3 = num_ave;
          hit_4 = num_ave2; 
        }
        if(file_flag == 1) 
          printf("frequency:%d %d %d %d %d %d %d %d %f %f %f || %d %d %d || %d %d\n",n_links[i],i*nstep,n_cells,i_getindex,j_getindex,R_len,ave5_1,ave5_2,rr,rm,rn,num_ave,num_ave2,max_locu,cutave_1,cutave_2);
     } 

     if((rate > 0.0025)&&(max_locu > 110000)&&(max_locu < (R_len - 110000)))
     {
       int idt = max_locu/nstep;
       int ave_left = 0;
       int ave_right = 0;
       int ave_edge = 0;

       for(j=0;j<idt;j++)
          ave_left = ave_left + n_links[j];
       for(j=idt;j<n_cells;j++)
          ave_right = ave_right + n_links[j]; 

       ave_left = ave_left/idt;
       ave_right = ave_right/(n_cells-idt);

       if(ave_left > ave_right)
       {
         if(ave_right == 0)
           ave_edge = 100;
         else
           ave_edge = ave_left/ave_right;
       }
       else
       {
         if(ave_left == 0)
           ave_edge = 100;
         else
           ave_edge = ave_right/ave_left;
       }
       d_score1 = max_rate;
       d_score2 = max_rate3;
       n_edges  = ave_edge;
       if((ave_edge >= 5)&&(max_rate > 0.6)&&(max_rate3 > 0.7))
       {
         printf("Break: %d %d %f %f %f %d || %d %d %f %f || %d %d\n",i_getindex,j_getindex,max_rate,max_rate2,max_rate3,ave_edge,R_len,max_locu,rate,all_rate,ave_left,ave_right);
         Break_Location(nSeq,R_len,max_locu,max_step,rd_locus,rd_index);
       } 
       printf("Max: %d %d %f %f %f %d || %d %d %f %f || %d %d\n",i_getindex,j_getindex,max_rate,max_rate2,max_rate3,ave_edge,R_len,max_locu,rate,all_rate,ave_left,ave_right);
     }
}


/* ========================================================================================= */
void Break_Location(int nSeq, int R_len, int m_loci, int m_step, int *rd_locus, int *rd_index)
/* ========================================================================================== */
{
     int i,j,k,stopflag,num_steps,num_hits,num_reads;
     int loci1,loci2,n_hits,m_index,m_locus;
     float rate,rate2,rate3,m_rates;

     loci1 = 0;
     loci2 = 0;

     m_index = 0;
     m_locus = 0;
     m_rates = 0.0;
     for(i=0;i<(nSeq-42);i++)
     {
        int ave_gap = 0;
        int ave_gap2 = 0;

        for(j=(i-40);j<=(i-1);j++)
        {
           if(j>=0)
             ave_gap = ave_gap+(rd_locus[j+1]-rd_locus[j]);
        }
        ave_gap = ave_gap/40;

        for(j=(i+2);j<(i+42);j++)
        {
           if(j<nSeq)
             ave_gap2 = ave_gap2+(rd_locus[j+1]-rd_locus[j]);
        }
        ave_gap2 = ave_gap2/40;

        if(ave_gap2 < ave_gap)
        {
          rate = ave_gap;
          if(ave_gap2 == 0)
            ave_gap2 = 1;
          rate = rate/ave_gap2; 
        }
        else
        {
          rate = ave_gap2;
          if(ave_gap == 0)
          {
            ave_gap = 1;
            rate = ave_gap2; 
          }
          else
            rate = rate/ave_gap; 
        }
        printf("Gap:%d-%d-%d %d %d %d || %d %d %f %d\n",i_getindex,j_getindex,i,rd_locus[i],ave_gap,rd_locus[i+1]-rd_locus[i],ave_gap,ave_gap2,rate,nSeq);
        if((rate > m_rates)&&(rd_locus[i] > 100000)&&(rd_locus[i] < (R_len-100000)))
        {
          m_index = i;
          m_locus = rd_locus[i];
          m_rates = rate;
        }
//       if(abs(loci1 - ave_gap) > 2000)
//          printf("SeeBreak:%d-%d-%d %d %d %d\n",i_getindex,nSeq,i,rd_locus[i],ave_gap,rd_locus[i+1]-rd_locus[i]);
          
        loci1 = ave_gap;
        loci2 = map_loci[i+1]-map_loci[i];
     }
     printf("SeeBreak:%d-%d-%d %d %f\n",i_getindex,nSeq,m_index,m_locus,m_rates);
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


