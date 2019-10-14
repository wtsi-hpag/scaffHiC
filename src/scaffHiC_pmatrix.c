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
static int *ctg_length,*hit_index;

/* SSAS default parameters   */
static int IMOD=0;
static int nContig=0;
static int file_flag = 1;
static int min_len = 100000;
static int max_ctg = 1000000;
static int mtg_length = 100000000;
static int n_depth = 60;
static int c_pairs = 100;
static float m_score = 200.0;
static float d_score = 0.8;
static float c_score = 1.0;
static float NL_score = 0.0;
static float NR_score = 0.0;
static float N_score = 0.25;
static float l_score = 1.0;
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
      printf("Usage: %s -depth 60 -score 200 <Input_HiC-alignment> <Output_matrix_agp>\n",argv[0]);

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
         sscanf(argv[++i],"%d",&min_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%f",&m_score);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-degree"))
       {
         sscanf(argv[++i],"%f",&l_score);
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
    if(system("ps aux | grep scaffHiC_pmatrix; date") == -1)
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
    while(fscanf(namef,"%s %s %d %d %s %d",R_Name[i],rdname,&hit_locus1[i],&hit_score[i],tempc1,&superlength[i])!=EOF)
    {
        st = rdname;
        ed = strrchr(rdname,'_');
        idt = atoi(ed+1);
        if(idt > max_ctg)
          max_ctg = idt;
        ctg_length[idt] = superlength[i];
        hit_index[i] = idt;
        i++;
    }
    fclose(namef);

    max_ctg = max_ctg + 1;
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
     int *ctg_list,*ptp_list,*ptn_list,*ctg_print;
     long num_cells,n_Bases;
     int num_hits,num_hit1,num_hit2,rcdex,rsize,rsize2,size_row,dlinks[5];
     int stopflag,offset,*ray,*dex;
     void ArraySort_Mix(int n, long *arr, int *brr);
     void ArraySort_float2(int n, float *arr, int *brr);
     char **DBname,*st,*ed;
     int **p_matrix,**PP_matrix,**PN_matrix,**PO_matrix,**s_matrix,**s2_matrix,**o_matrix,**r_matrix,**rc0_matrix,**rc1_matrix,**rc2_matrix,**rc3_matrix,**rc_matrix,**rcdex_matrix;
     float rate,c_score1,c_score2,nl_score,nr_score,*Dis_index,*Dis_ratia1,*Dis_ratia2,*Dis_score1,*Dis_score2,*Dnd_score1,*Dnd_score2,*DD_score1,*DD_score2;
     int **imatrix(long nrl,long nrh,long ncl,long nch);
     float **fmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     void ArraySort_Int2(int n, int *arr, int *brr);
     void Reads_Distribution(int n, int m, int *arr, int *brr);
     void Normal_Distribution(int n, int m, int *arr, int *brr);
     void Direction_Ratio(int n,int *arr); 
     int *ctg_score1,*ctg_score2,*ctg_mapp1,*ctg_mapp2,*ctg_join1,*ctg_join2,*ctg_idex1,*ctg_idex2,*ctg_mask,*ctg_rcdex1;
     int *p_index,*p_rcdex,*p_masks,*p_score,*p_lists;
     int *ctg_rcoo1,*ctg_rcoo2,*ctg_part1,*ctg_part2,*ctg_patnum;
     int *ctg_output,*ctg_hitnum,*ctg_used,*ctg_links,*ctg_oodex1,*ctg_oodex2,*ctg_rcindex,*ctg_mpindex,*ctg_outrc;
     int *hit_linksdex,*link_locus,*link_locu2,*link_locu3,*head_locus;
     int n_length = min_len;

     rsize = max_ctg+10; 
     n_blocks = rsize*rsize;
     len_thresd = 25000000;
 
     if((ctg_idex1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_idex1\n");
       exit(1);
     }
     if((ctg_idex2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_idex2\n");
       exit(1);
     }
     if((ctg_rcdex1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_mapp1\n");
       exit(1);
     }
     if((ctg_mapp1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_mapp1\n");
       exit(1);
     }
     if((ctg_mapp2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_mapp2\n");
       exit(1);
     }
     if((ctg_part1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_part1\n");
       exit(1);
     }
     if((ctg_part2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_part2\n");
       exit(1);
     }
     if((ctg_score1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score1\n");
       exit(1);
     }
     if((ctg_score2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }
     if((ctg_join1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_join1\n");
       exit(1);
     }
     if((ctg_join2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_join2\n");
       exit(1);
     }
     if((ctg_mask = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_mask\n");
       exit(1);
     }
     if((ctg_oodex1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_oodex1\n");
       exit(1);
     }
     if((ctg_oodex2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_oodex2\n");
       exit(1);
     }
     if((ctg_rcoo1 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_rcoo1\n");
       exit(1);
     }
     if((ctg_rcoo2 = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_rcoo2\n");
       exit(1);
     }
     if((ctg_hitnum = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_hitnum\n");
       exit(1);
     }
     if((ctg_patnum = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_pattnum\n");
       exit(1);
     }
     if((ctg_output = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_output\n");
       exit(1);
     }
     if((ctg_used = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_used\n");
       exit(1);
     }
     if((ctg_links = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_links\n");
       exit(1);
     }
     if((ctg_rcindex = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_rcindex\n");
       exit(1);
     }
     if((ctg_outrc = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_outrc\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((ptp_list = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ptp_list\n");
       exit(1);
     }
     if((ptn_list = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ptn_list\n");
       exit(1);
     }
     if((ctg_print = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_print\n");
       exit(1);
     }
     if((ctg_mpindex = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_mpindex\n");
       exit(1);
     }
     if((Dis_index = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_index\n");
       exit(1);
     }
     if((Dis_ratia1 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_index\n");
       exit(1);
     }
     if((Dis_ratia2 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_index\n");
       exit(1);
     }
     if((DD_score1 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_score1\n");
       exit(1);
     }
     if((DD_score2 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_score1\n");
       exit(1);
     }
     if((Dis_score1 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_score1\n");
       exit(1);
     }
     if((Dis_score2 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dis_score2\n");
       exit(1);
     }
     if((Dnd_score1 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dnd_score1\n");
       exit(1);
     }
     if((Dnd_score2 = (float *)calloc(rsize,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - Dnd_score2\n");
       exit(1);
     }

     if((p_index = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - p_index\n");
       exit(1);
     }
     if((p_rcdex = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - p_rcdex\n");
       exit(1);
     }
     if((p_score = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - p_score\n");
       exit(1);
     }
     if((p_masks = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - p_mask\n");
       exit(1);
     }
     if((p_lists = (int *)calloc(rsize,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - p_lists\n");
       exit(1);
     }

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
     p_matrix=imatrix(0,rsize,0,rsize);
     PP_matrix=imatrix(0,rsize,0,100);
     PN_matrix=imatrix(0,rsize,0,100);
     PO_matrix=imatrix(0,rsize,0,100);
     s_matrix=imatrix(0,rsize,0,rsize);
     s2_matrix=imatrix(0,rsize,0,rsize);
     o_matrix=imatrix(0,rsize,0,rsize);
     r_matrix=imatrix(0,rsize,0,rsize);
     rc0_matrix=imatrix(0,rsize,0,rsize);
     rc1_matrix=imatrix(0,rsize,0,rsize);
     rc2_matrix=imatrix(0,rsize,0,rsize);
     rc3_matrix=imatrix(0,rsize,0,rsize);
     rc_matrix=imatrix(0,rsize,0,rsize);
     rcdex_matrix=imatrix(0,rsize,0,rsize);
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
     for(i=0;i<max_ctg;i++)
     {
        ctg_part1[i] = -1;
        ctg_part2[i] = -1;
        ctg_mapp1[i] = -1;
        ctg_mapp2[i] = -1;
        ctg_rcdex1[i] = -1;
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

     n_Bases = 0;
     for(i=0;i<max_ctg;i++)
        n_Bases = n_Bases + ctg_length[i];

     printf("Total number of contigs: %d\n",max_ctg);
     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
           hit_rddex[j] = j;
        Dis_index[i] = 0.0;
        ArraySort2_Int2(max_ctg,r_matrix[i],hit_rddex);
 
        memset(ctg_rcindex,0,4*max_ctg); 

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

           ctg_rcindex[j] = rcindex; 

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
           c_score1 = 1.0;
           if((rat1 > 1.8)||(rat2 > 1.8)||(ctg_length[idi] < 1000000))
           {
             i_getindex = idi;
             j_getindex = idj;  
             Reads_Distribution(n_pairs,ctg_length[idi],ray+offset1,dex+offset1);
             printf("Distribute: %d %d %d %d %f %f || %d %d\n",halflen_1,halflen_2,half_len1,half_len2,d_score,c_score,n_pairs,c_pairs);
             c_score1 = c_score;
             c1_pairs = c_pairs;
           }
           c_score2 = 1.0;
//           if((rat1 > 2.0)||(rat2 > 2.0))
           if((rat1 > 1.8)||(rat2 > 1.8)||(ctg_length[idj] < 1000000))
           {
             i_getindex = idj;
             j_getindex = idi;  
             Reads_Distribution(n_pairs,ctg_length[idj],ray+offset2,dex+offset2);
             printf("Distribute: %d %d %d %d %f %f || %d %d\n",halflen_1,halflen_2,half_len1,half_len2,d_score,c_score,n_pairs,c_pairs);
             c_score2 = c_score;
             c2_pairs = c_pairs;
           }

           if(c1_pairs > c2_pairs)
             c_pairs = c2_pairs;
           else
             c_pairs = c1_pairs; 

           c_pairs = n_pairs;
           rr1 = n_Bases;
           rr1 = rr1/nSeq;
           rr1 = rr1*c_pairs*1000.0; 
           rr1 = rr1/half_len1;
           rr1 = rr1*(rat1 - 1.0);
//           rr1 = rr1/half_len1;

           rr2 = n_Bases;
           rr2 = rr2/nSeq;
           rr2 = rr2*c_pairs*1000.0; 
           rr2 = rr2/half_len2;
           rr2 = rr2*(rat2 - 1.0);
//           rr2 = rr2/half_len2;

           nl_score = 0.0;
           nr_score = 0.0;
           if(ctg_length[idi] > ctg_length[idj])
           {
             rat_size = ctg_length[idi];
             rat_size = rat_size/ctg_length[idj];
             if(rat_size > 30.0)
             {
               i_getindex = idi;
               j_getindex = idj;  
               Normal_Distribution(n_pairs,ctg_length[idi],ray+offset1,dex+offset1);
               nl_score = NL_score;
               nr_score = NR_score;
//               printf("Distribute: %d %d %d %d | %f %f\n",n_pairs,idi,idj,ctg_length[idi],nl_score,nr_score);
             } 
           }
           else
           {
             rat_size = ctg_length[idj];
             rat_size = rat_size/ctg_length[idi];
             if(rat_size > 30.0)
             {
               i_getindex = idi;
               j_getindex = idj;  
               Normal_Distribution(n_pairs,ctg_length[idj],ray+offset2,dex+offset2);
               nl_score = NL_score;
               nr_score = NR_score;
//               printf("Distribute: %d %d %d %d | %f %f\n",n_pairs,idi,idj,ctg_length[idi],nl_score,nr_score);
             } 
           }

           DD_score1[j] = rr1;
           DD_score2[j] = rr2;
           Dis_score1[j] = c_score1;
           Dis_score2[j] = c_score2;
           Dnd_score1[j] = nl_score;
           Dnd_score2[j] = nr_score;
           Dis_ratia1[j] = rat1;
           Dis_ratia2[j] = rat2;

           Dis_index[j] = rr1*rr2; 
           printf("%6d | %6d %d %d | %d %d %f %f %f %f | %f %f %f %4d %4d %4d %4d | %f %f\n",ctg_length[idd],r_matrix[i][j],hit_rddex[j],rcindex,halflen_1,halflen_2,rat1,rat2,c_score1,c_score2,rr1,rr2,rr1*rr2,rc0_matrix[i][idd],rc1_matrix[i][idd],rc2_matrix[i][idd],rc3_matrix[i][idd],nl_score,nr_score);
        }

        if((ctg_length[i]>=n_length))
        {
          int OO_index1 = 0;
          int OO_index2 = 0;
          int hitmax1 = 0;
          int hitmax2 = 0;
          int ctgmax1 = 0;
          int ctgmax2 = 0;
          int idi,idj,disidd;
          float disdex = 0.0;
          float disdex2 = 0.0;
          float rr,rr2,rr3,M_score;
         
            OO_index1 = 0;
            for(k=0;k<max_ctg;k++)
            {
               float DS = 0.0;
               if(DD_score1[k] >= DD_score2[k])
               {
                 DS = DD_score1[k];
                 if(DD_score2[k] == 0.0)
                   DS = 1000;
                 else
                   DS = DS/DD_score2[k];
               }
               else
               {
                 DS = Dis_score2[k];
                 if(DD_score1[k] == 0.0)
                   DS = 1000;
                 else
                   DS = DS/DD_score1[k];
               }

//               if(i==41)
//                 printf("Dis: %d %f %f %f\n",k,DD_score1[k],DD_score2[k],DS);
               if((ctg_length[hit_rddex[k]]>=n_length)&&(disdex < Dis_index[k])&&(ctg_rcindex[k]<2)&&(Dis_ratia1[k] > 1.02)&&(Dis_ratia2[k] > 1.02)&&(Dis_score1[k] >= d_score)&&(Dis_score2[k] >= d_score)&&(Dnd_score1[k]<=N_score)&&(Dnd_score2[k]<=N_score)&&(DS < 50))
               {
                 disdex = Dis_index[k];
                 OO_index1 = ctg_rcindex[k];
                 disidd = k;
               } 
            }
           
            if((ctg_length[i] < len_thresd)&&(ctg_length[hit_rddex[disidd]]< len_thresd))
            {
              M_score = m_score;
            }
            else
            {
              float set1 = ctg_length[i];
              float set2 = ctg_length[hit_rddex[disidd]];
              set1 = set1/len_thresd;
              set2 = set2/len_thresd;
              M_score = (set1+set2)*m_score*l_score;
            if(i==5||i==254)
              printf("score00: %d %f %f %f %f\n",i,M_score,disdex,set1,set2); 
            }

            if(i==5||i==254)
              printf("score1: %d %f %f\n",i,M_score,disdex);
            dlinks[0] = rc0_matrix[i][hit_rddex[disidd]];
            dlinks[1] = rc1_matrix[i][hit_rddex[disidd]];
            dlinks[2] = rc2_matrix[i][hit_rddex[disidd]];
            dlinks[3] = rc3_matrix[i][hit_rddex[disidd]];
            dlinks[4] = 0;
            Direction_Ratio(i,dlinks); 
            if((disdex > M_score)&&(dlinks[4] == 1)) 
            {
               if(dlinks[0] == rc0_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 1;
                 rcdex_matrix[hit_rddex[disidd]][i] = 1;
               }
               else if(dlinks[0] == rc1_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 2;
                 rcdex_matrix[hit_rddex[disidd]][i] = 3;
               } 
               else if(dlinks[0] == rc2_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 3;
                 rcdex_matrix[hit_rddex[disidd]][i] = 2;
               } 
               else if(dlinks[0] == rc3_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 4;
                 rcdex_matrix[hit_rddex[disidd]][i] = 4;
               } 
               if((ctg_length[i] >= mtg_length)&&(ctg_length[hit_rddex[disidd]] >= mtg_length))
                 printf("mapp: 1 max %d %d %d || %f %d %d %f %f || %d %d %d | %f %f\n",i,hit_rddex[disidd],r_matrix[i][disidd],disdex,OO_index1,ctg_oodex1[i],Dis_ratia1[disidd],Dis_ratia2[disidd],ctg_length[i],ctg_length[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],Dnd_score1[disidd],M_score);
               else
                 printf("mapp: 1 %d %d %d || %f %d %d %f %f || %d %d %d | %f %f\n",i,hit_rddex[disidd],r_matrix[i][disidd],disdex,OO_index1,ctg_oodex1[i],Dis_ratia1[disidd],Dis_ratia2[disidd],ctg_length[i],ctg_length[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],Dnd_score1[disidd],M_score);
               ctg_mapp1[i] = hit_rddex[disidd];
               ctg_score1[i] = disdex;
               ctg_rcoo1[i] = OO_index1;
/*         add new lines   */
               if((rcdex_matrix[i][hit_rddex[disidd]] == 1)||(rcdex_matrix[i][hit_rddex[disidd]] == 4))
                 rc_matrix[i][hit_rddex[disidd]] = 1;
               else
                 rc_matrix[i][hit_rddex[disidd]] = 0;
               rc_matrix[hit_rddex[disidd]][i] = rc_matrix[i][hit_rddex[disidd]];
 
               PP_matrix[i][ptp_list[i]] = hit_rddex[disidd];
               ptp_list[i]++;
               PN_matrix[hit_rddex[disidd]][ptn_list[hit_rddex[disidd]]] = i;
               if((hit_rddex[disidd] == 257)||(i == 257))
                 printf("order: 1 %d %d %d %d\n",i,ptn_list[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],rcdex_matrix[hit_rddex[disidd]][i]);
               PO_matrix[hit_rddex[disidd]][ptn_list[hit_rddex[disidd]]] = rcdex_matrix[hit_rddex[disidd]][i];
               ptn_list[hit_rddex[disidd]]++;
               ctg_hitnum[i]++;
            }


            OO_index2 = 0;
            for(k=0;k<max_ctg;k++)
            {
               float DS = 0.0;
               if(DD_score1[k] >= DD_score2[k])
               {
                 DS = DD_score1[k];
                 if(DD_score2[k] == 0.0)
                   DS = 1000;
                 else
                   DS = DS/DD_score2[k];
               }
               else
               {
                 DS = DD_score2[k];
                 if(DD_score1[k] == 0.0)
                   DS = 1000;
                 else
                   DS = DS/DD_score1[k];
               }
               if((ctg_length[hit_rddex[k]]>=n_length)&&(disdex2 < Dis_index[k])&&(ctg_rcindex[k]>=2)&&(Dis_ratia1[k] > 1.02)&&(Dis_ratia2[k] > 1.02)&&(Dis_score1[k] >= d_score)&&(Dis_score2[k] >= d_score)&&(Dnd_score1[k]<=N_score)&&(Dnd_score2[k]<=N_score)&&(DS < 50))
               {
                 disdex2 = Dis_index[k];
                 OO_index2 = ctg_rcindex[k];
                 disidd = k;
               } 
            }
            if((ctg_length[i] < len_thresd)&&(ctg_length[hit_rddex[disidd]]< len_thresd))
            {
              M_score = m_score;
            }
            else
            {
              float set1 = ctg_length[i];
              float set2 = ctg_length[hit_rddex[disidd]];
              set1 = set1/len_thresd;
              set2 = set2/len_thresd;
              M_score = (set1+set2)*m_score*l_score;
            } 
            dlinks[0] = rc0_matrix[i][hit_rddex[disidd]];
            dlinks[1] = rc1_matrix[i][hit_rddex[disidd]];
            dlinks[2] = rc2_matrix[i][hit_rddex[disidd]];
            dlinks[3] = rc3_matrix[i][hit_rddex[disidd]];
            dlinks[4] = 0;
            Direction_Ratio(i,dlinks); 
            if((disdex2 > M_score)&&(dlinks[4] == 1)) 
            {
               if(dlinks[0] == rc0_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 1;
                 rcdex_matrix[hit_rddex[disidd]][i] = 1;
               }
               else if(dlinks[0] == rc1_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 2;
                 rcdex_matrix[hit_rddex[disidd]][i] = 3;
               } 
               else if(dlinks[0] == rc2_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 3;
                 rcdex_matrix[hit_rddex[disidd]][i] = 2;
               } 
               else if(dlinks[0] == rc3_matrix[i][hit_rddex[disidd]])
               {
                 rcdex_matrix[i][hit_rddex[disidd]] = 4;
                 rcdex_matrix[hit_rddex[disidd]][i] = 4;
               } 
               if((ctg_length[i] >= mtg_length)&&(ctg_length[hit_rddex[disidd]] >= mtg_length))
                 printf("mapp: 2 max %d %d %d || %f %d %d %f %f || %d %d %d | %f %f\n",i,hit_rddex[disidd],r_matrix[i][disidd],disdex2,OO_index2,ctg_oodex2[i],Dis_ratia1[disidd],Dis_ratia2[disidd],ctg_length[i],ctg_length[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],Dnd_score1[disidd],M_score);
               else
                 printf("mapp: 2 %d %d %d || %f %d %d %f %f || %d %d %d | %f %f\n",i,hit_rddex[disidd],c_pairs,disdex2,OO_index2,ctg_oodex2[i],Dis_ratia1[disidd],Dis_ratia2[disidd],ctg_length[i],ctg_length[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],Dnd_score1[disidd],M_score);
//                 printf("mapp: 2 %d %d %d || %f %d %d %f %f || %d %d %d | %f %f\n",i,hit_rddex[disidd],r_matrix[i][disidd],disdex2,OO_index2,ctg_oodex2[i],Dis_ratia1[disidd],Dis_ratia2[disidd],ctg_length[i],ctg_length[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],Dnd_score1[disidd],M_score);
               ctg_mapp2[i] = hit_rddex[disidd];
               ctg_score2[i] = disdex2;
               ctg_rcoo2[i] = OO_index2;
               ctg_rcdex1[i] = OO_index2;

/*         add new lines   */
               if((rcdex_matrix[i][hit_rddex[disidd]] == 1)||(rcdex_matrix[i][hit_rddex[disidd]] == 4))
                 rc_matrix[i][hit_rddex[disidd]] = 1;
               else
                 rc_matrix[i][hit_rddex[disidd]] = 0;
               rc_matrix[hit_rddex[disidd]][i] = rc_matrix[i][hit_rddex[disidd]];
 
               PP_matrix[i][ptp_list[i]] = hit_rddex[disidd];
               ptp_list[i]++;
               PN_matrix[hit_rddex[disidd]][ptn_list[hit_rddex[disidd]]] = i;
               if((hit_rddex[disidd] == 257)||(i == 257))
                 printf("order: 2 %d %d %d %d\n",i,ptn_list[hit_rddex[disidd]],rcdex_matrix[i][hit_rddex[disidd]],rcdex_matrix[hit_rddex[disidd]][i]);
               PO_matrix[hit_rddex[disidd]][ptn_list[hit_rddex[disidd]]] = rcdex_matrix[hit_rddex[disidd]][i];
               ptn_list[hit_rddex[disidd]]++;
               ctg_hitnum[i]++;
            }
        } 
     }

     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
        {
           if((rcdex_matrix[j][i] > 0)&&(rcdex_matrix[i][j] == 0))
           {
             if(rcdex_matrix[j][i] == 1)
               rcdex_matrix[i][j] = 4;
             else if(rcdex_matrix[j][i] == 2)
               rcdex_matrix[i][j] = 3; 
             else if(rcdex_matrix[j][i] == 3)
               rcdex_matrix[i][j] = 2; 
             else if(rcdex_matrix[j][i] == 4)
               rcdex_matrix[i][j] = 1; 
           }
           if((rcdex_matrix[j][i] == 0)&&(rcdex_matrix[i][j] > 0))
           {
             if(rcdex_matrix[i][j] == 1)
               rcdex_matrix[j][i] = 4;
             else if(rcdex_matrix[i][j] == 2)
               rcdex_matrix[j][i] = 3; 
             else if(rcdex_matrix[i][j] == 3)
               rcdex_matrix[j][i] = 2; 
             else if(rcdex_matrix[i][j] == 4)
               rcdex_matrix[j][i] = 1; 
           }
        }
     }

     for(k=0;k<max_ctg;k++)
     {
        if(ctg_hitnum[k] == 1)
        {
          int idi = ctg_mapp1[k];
          if(ctg_mapp1[k] >= 0)
          {
            p_matrix[k][ctg_list[k]] = ctg_mapp1[k];
            s_matrix[k][ctg_list[k]] = ctg_score1[k];
            s2_matrix[k][ctg_mapp1[k]] = ctg_score1[k];
            s2_matrix[ctg_mapp1[k]][k] = ctg_score1[k];
            o_matrix[k][ctg_list[k]] = ctg_rcoo1[k];
            ctg_list[k]++; 
          }
          else if(ctg_mapp2[k] >= 0)
          {
            p_matrix[k][ctg_list[k]] = ctg_mapp2[k];
            s_matrix[k][ctg_list[k]] = ctg_score2[k];
            s2_matrix[k][ctg_mapp2[k]] = ctg_score2[k];
            s2_matrix[ctg_mapp2[k]][k] = ctg_score2[k];
            o_matrix[k][ctg_list[k]] = ctg_rcoo2[k];
            ctg_list[k]++;
          } 
          for(j=0;j<max_ctg;j++)
          {
             if((ctg_mapp1[j] == k)||(ctg_mapp2[j] == k))
             {
               int idi = ctg_hitnum[ctg_mapp1[j]];
               int idj = ctg_hitnum[ctg_mapp2[j]];
               printf("find the missing link: %d %d %d %d %d %d\n",k,j,ctg_mapp1[j],ctg_mapp2[j],idi,idj);
             }
          }
        }
        else if(ctg_hitnum[k] == 2)
        {
          p_matrix[k][ctg_list[k]] = ctg_mapp1[k];
          s_matrix[k][ctg_list[k]] = ctg_score1[k];
          s2_matrix[k][ctg_mapp1[k]] = ctg_score1[k];
          s2_matrix[ctg_mapp1[k]][k] = ctg_score1[k];
          o_matrix[k][ctg_list[k]] = ctg_rcoo1[k];
          ctg_list[k]++;
          p_matrix[k][ctg_list[k]] = ctg_mapp2[k];
          s_matrix[k][ctg_list[k]] = ctg_score2[k];
          s2_matrix[k][ctg_mapp2[k]] = ctg_score2[k];
          s2_matrix[ctg_mapp2[k]][k] = ctg_score2[k];
          o_matrix[k][ctg_list[k]] = ctg_rcoo2[k];
          ctg_list[k]++;
        }         
     }

     
     for(k=0;k<max_ctg;k++)
     {
        int num_partner = 0;
        int p_idi = -1;
        int p_idk = -1;
        int pscore1 = 0;         
        int pscore2 = 0;         
        int plists1 = 0;
        int plists2 = 0;

        for(j=0;j<max_ctg;j++)
        {
           for(i=0;i<ctg_list[j];i++)
           {
              if(p_matrix[j][i] == k)
              {
                p_index[num_partner] = j;
                p_rcdex[num_partner] = o_matrix[j][i];
                p_score[num_partner] = s_matrix[j][i];
                p_lists[num_partner] = i;
                num_partner++;
              }
           }

           for(i=0;i<ctg_list[j];i++)
           {
              if(p_matrix[j][i] == k)
                 printf("Missing link: %d %d | %d %d %d %d | %d %d || %d %d %d %d || %d %d\n",k,j,i,p_matrix[j][i],s_matrix[j][i],o_matrix[j][i],ctg_length[k],ctg_length[j],rc0_matrix[k][j],rc1_matrix[k][j],rc2_matrix[k][j],rc3_matrix[k][j],ctg_mapp1[k],ctg_mapp2[k]);
           }
        }

        if(num_partner > 0)
        {
          for(i=0;i<num_partner;i++)
          {
             if((p_rcdex[i]%2) == 0)
             {
               if(p_score[i] > pscore1)
               {
                 p_idi = i;
                 pscore1 = p_score[i];
                 plists1 = p_lists[i];
               } 
             }
             else
             {
               if(p_score[i] > pscore2)
               {
                 p_idk = i;
                 pscore2 = p_score[i];
                 plists2 = p_lists[i];
               }
             } 
          }  
        }

        if(p_idi >= 0)
        {
          j = p_index[p_idi];
          i = plists1;
          if((o_matrix[j][i] == 1)||(o_matrix[j][i] == 2))
            ctg_oodex1[k] = 0;
          else
            ctg_oodex1[k] = 1;
          rc_matrix[k][j] = ctg_oodex1[k];
          rc_matrix[j][k] = ctg_oodex1[k];
          ctg_part1[k] = j;
          ctg_patnum[k]++; 
          ctg_links[k]++; 
          printf("Partner_1: %d %d %d | %d %d | %d %d || %d %d %d %d %d\n",k,j,i,s_matrix[j][i],o_matrix[j][i],ctg_length[k],ctg_length[j],rc0_matrix[k][i],rc1_matrix[k][i],rc2_matrix[k][i],rc3_matrix[k][i],ctg_mapp1[k]);
        }
        if(p_idk >= 0)
        {
          j = p_index[p_idk];
          i = plists2;
          if((o_matrix[j][i] == 1)||(o_matrix[j][i] == 2))
            ctg_oodex2[k] = 0;
          else
            ctg_oodex2[k] = 1;
          rc_matrix[k][j] = ctg_oodex2[k];
          rc_matrix[j][k] = ctg_oodex2[k];
          ctg_part2[k] = j;
          ctg_patnum[k]++; 
          ctg_links[k]++; 
          if(j != ctg_mapp2[k])
          {
            ctg_part1[j] = ctg_mapp2[k]; 
                 printf("PPPpartner: %d %d | %d %d\n",k,j,ctg_mapp1[k],ctg_mapp2[k]);
          }
          printf("Partner_2: %d %d %d | %d %d | %d %d || %d %d %d %d %d\n",k,j,i,s_matrix[j][i],o_matrix[j][i],ctg_length[k],ctg_length[j],rc0_matrix[k][j],rc1_matrix[k][j],rc2_matrix[k][j],rc3_matrix[k][j],ctg_mapp2[k]);
        }
        for(i=0;i<ctg_list[k];i++)
        {
              printf("link: %d %d %d %d %d\n",k,i,ctg_list[k],p_matrix[k][i],s_matrix[k][i]);
        }
     }


     for(k=0;k<max_ctg;k++)
     {
        for(i=0;i<ctg_list[i];i++)
           printf("k-428: %d %d %d %d | %d %d | %d %d\n",k,i,ctg_list[k],p_matrix[k][i],ctg_mapp1[k],ctg_mapp2[k],ctg_part1[k],ctg_part2[k]);
     }
     for(k=0;k<max_ctg;k++)
     {
        if(ptp_list[k] == 2)
        {
          if(ptn_list[k] == 0)
          {
            ctg_print[k] = 1;
            ctg_print[PP_matrix[k][0]] = 1;
            ctg_print[PP_matrix[k][1]] = 1;
            printf("next: %d %d %d %d %d\n",k,ctg_part1[k],ctg_part2[k],PP_matrix[k][0],PP_matrix[k][1]);
            ctg_part1[k] = PP_matrix[k][0];
            ctg_part2[k] = PP_matrix[k][1];
          }
          else if(ptn_list[k] == 1)
          {
            if((PN_matrix[k][0] == PP_matrix[k][0])||(PN_matrix[k][0] == PP_matrix[k][1]))
            {
              ctg_print[k] = 1;
              ctg_print[PP_matrix[k][0]] = 1;
              ctg_print[PP_matrix[k][1]] = 1;
            printf("next: %d %d %d %d %d\n",k,ctg_part1[k],ctg_part2[k],PP_matrix[k][0],PP_matrix[k][1]);
              ctg_part1[k] = PP_matrix[k][0];
              ctg_part2[k] = PP_matrix[k][1];
            }
          }
          else if(ptn_list[k] == 2)
          {
            if((PN_matrix[k][0] == PP_matrix[k][0])&&(PN_matrix[k][1] == PP_matrix[k][1]))
            {
              ctg_print[k] = 1;
              ctg_print[PP_matrix[k][0]] = 1;
              ctg_print[PP_matrix[k][1]] = 1;
            printf("next: %d %d %d %d %d\n",k,ctg_part1[k],ctg_part2[k],PP_matrix[k][0],PP_matrix[k][1]);
              ctg_part1[k] = PP_matrix[k][0];
              ctg_part2[k] = PP_matrix[k][1];
            }
            else if((PN_matrix[k][0] == PP_matrix[k][1])&&(PN_matrix[k][1] == PP_matrix[k][0]))
            {
              ctg_print[k] = 1;
              ctg_print[PP_matrix[k][0]] = 1;
              ctg_print[PP_matrix[k][1]] = 1;
            printf("next: %d %d %d %d %d\n",k,ctg_part1[k],ctg_part2[k],PP_matrix[k][0],PP_matrix[k][1]);
              ctg_part1[k] = PP_matrix[k][0];
              ctg_part2[k] = PP_matrix[k][1];
            }
          }
        }
        else if(ptp_list[k] == 1)
        {
          if(ptn_list[k] == 1)
          {
            if((ptn_list[k] == 1)&&(PP_matrix[k][0]==PN_matrix[k][0]))
            {
              ctg_print[k] = 1;
              ctg_print[PP_matrix[k][0]] = 1;
            printf("next: %d %d %d %d %d\n",k,ctg_part1[k],ctg_part2[k],PP_matrix[k][0],100);
              ctg_part1[k] = PP_matrix[k][0];
            } 
          }
        }
     }

     for(k=0;k<max_ctg;k++)
     {
        printf("PPmatrix: %d %d | %d %d | %d %d %d || ",k,ptp_list[k],ctg_part1[k],ctg_part2[k],ptn_list[k],ctg_patnum[k],ctg_print[k]);
//        printf("PPmatrix: %d %d | %d %d %d || ",k,ptp_list[k],ptn_list[k],ctg_patnum[k],ctg_print[k]);
        for(i=0;i<ptp_list[k];i++)
           printf("%d ",PP_matrix[k][i]);
        printf("|| ");
        for(i=0;i<ptn_list[k];i++)
           printf("%d ",PN_matrix[k][i]);
        printf("& ");
        for(i=0;i<ptn_list[k];i++)
           printf("%d ",PO_matrix[k][i]);
        printf("\n");

     }

     for(k=0;k<max_ctg;k++)
     {
        if((ptn_list[k] == 1)&&(ptp_list[k] == 2))
        {
          if(PN_matrix[k][0] == PP_matrix[k][0])
          {
            ctg_print[k] = 1;
            ctg_print[PP_matrix[k][1]] = 1;
            PN_matrix[k][1] = PP_matrix[k][1];
            PO_matrix[k][1] = rcdex_matrix[k][PP_matrix[k][1]];
            ptn_list[k]++;
          }
          else if(PN_matrix[k][0] == PP_matrix[k][1])
          {
            ctg_print[k] = 1;
            ctg_print[PP_matrix[k][0]] = 1;
            PN_matrix[k][1] = PP_matrix[k][0];
            PO_matrix[k][1] = rcdex_matrix[k][PP_matrix[k][0]];
            ptn_list[k]++;
          }
        }
        else if((ptn_list[k] == 0)&&(ptp_list[k] == 1))
        {
            ctg_print[k] = 1;
            ctg_print[PP_matrix[k][0]] = 1;
            PN_matrix[k][0] = PP_matrix[k][0];
            PO_matrix[k][0] = rcdex_matrix[k][PP_matrix[k][0]];
            ptn_list[k]++;
        }
        else if((ptn_list[k] == 0)&&(ptp_list[k] == 2))
        {
            ctg_print[k] = 1;
            ctg_print[PP_matrix[k][0]] = 1;
            ctg_print[PP_matrix[k][1]] = 1;
            PN_matrix[k][0] = PP_matrix[k][0];
            PO_matrix[k][0] = rcdex_matrix[k][PP_matrix[k][0]];
            ptn_list[k]++;
            PN_matrix[k][1] = PP_matrix[k][1];
            PO_matrix[k][1] = rcdex_matrix[k][PP_matrix[k][1]];
            ptn_list[k]++;
        }
     }

     for(k=0;k<max_ctg;k++)
     {
        printf("PPmatrix2: %d %d | %d %d | %d %d %d || ",k,ptp_list[k],ctg_part1[k],ctg_part2[k],ptn_list[k],ctg_patnum[k],ctg_print[k]);
//        printf("PPmatrix: %d %d | %d %d %d || ",k,ptp_list[k],ptn_list[k],ctg_patnum[k],ctg_print[k]);
        for(i=0;i<ptp_list[k];i++)
           printf("%d ",PP_matrix[k][i]);
        printf("|| ");
        for(i=0;i<ptn_list[k];i++)
           printf("%d ",PN_matrix[k][i]);
        printf("& ");
        for(i=0;i<ptn_list[k];i++)
           printf("%d ",PO_matrix[k][i]);
        printf("\n");

     }

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

          n_scaff = 0;
          for(k=0;k<max_ctg;k++)
          {
            int tbase = 0;
            int print_tag = 0;

            if(ptp_list[k] == 1)
            {
              if(ptn_list[k] == 1)
              {
                if((ptn_list[k] == 1)&&(PP_matrix[k][0]==PN_matrix[k][0]))
                  print_tag = 1;
              }
            }
//             if(k==i_getindex)
               printf("===========\n");
               printf("hhh: %d %d %d %d\n",k,ctg_patnum[k],ctg_output[k],ptn_list[k]);
            if((ptn_list[k]==0)&&(ctg_output[k]==0))
            {
              printf("supercontig1: tarseq_%d %d %d\n",k,k,ctg_patnum[k]);
              fprintf(namef,"contig-1:%08d %d %d %d 0 0\n",0,n_scaff,k,ctg_length[k]);
              printf("contig-1: %d %d %d 0 0\n",n_scaff,k,ctg_length[k]);
              tbase = ctg_length[k];
              printf("bases: %d %d %d\n",k,n_scaff,tbase);
              ctg_output[k] = 1;
              n_scaff++;
            }
            else if((print_tag==1)&&(ptn_list[k]==1)&&(ctg_output[k]==0))
            {
             int idd = k;
             int idk = k;
             int stopflag=0;
             int num_loops=0;

//             if(k==i_getindex)
//             while(((ctg_part1[idd] >= 0)||(ctg_part2[idd] >= 0))&&((stopflag == 0)&&(ctg_links[idk]>0)&&(ctg_mask[idk]==0)&&(idk >= 0)&&(idd >= 0)))
//             while((ctg_print[idd] == 1)&&((stopflag == 0)&&(ctg_links[idk]>0)&&(ctg_mask[idk]==0)&&(idk >= 0)&&(idd >= 0)))
             while((ctg_print[idd] == 1)&&((stopflag == 0)&&(ctg_mask[idk]==0)&&(idk >= 0)&&(idd >= 0)))
             {
               int rc_idk = -1;
               int rc_idi = 0;
               printf("xxx: %d %d %d %d %d %d %d %d\n",k,idd,idk,ptn_list[k],ctg_used[k],ctg_used[idd],ctg_part1[idd],ctg_part2[idd]);
               if(ctg_used[idk] == 0)
               {
             printf("hits: %d %d %d %d %d %d %d %d\n",k,idd,idk,ptn_list[idd],ctg_used[k],ctg_used[idd],ctg_part1[idd],ctg_part2[idd]);
                 if(ctg_used[idd] == 0)
                 {
                   int min_len = 2000000000;
                   int min_idd = 0;
                   int cono = 0;

                   if(ptn_list[idd] == 1)
                   {
                     idk = PN_matrix[idd][0];
                     if((PO_matrix[idd][0] >= 2)&&(PO_matrix[idd][0] <= 3))
                       rc_idk = 0;
                     else
                       rc_idk = 1;
                   }
                   else
                   {
                     for(m=0;m<ptn_list[idd];m++)
                     {
                        int rdex = PN_matrix[idd][m];

                        if(ctg_used[rdex] == 1)
                        {
                          if((PO_matrix[idd][m]%2) == 0)
                            cono = 2;
                          else
                            cono = 1; 
                        }
             printf("hits1: %d %d %d %d %d %d %d %d\n",m,idd,idk,ptn_list[idd],ctg_used[rdex],PO_matrix[idd][m],rdex,cono);
                     }
             printf("hits2: %d %d %d %d %d %d %d %d\n",k,idd,idk,ptn_list[idd],ctg_used[k],ctg_used[idd],ctg_part1[idd],cono);
                     for(m=0;m<ptn_list[idd];m++)
                     {
                        int rdex = PN_matrix[idd][m];
                        int cdex = ctg_part2[rdex];
                        int mono;
                      
                        mono = 0;
                      printf("hhh-: %d %d %d %d %d | %d %d || %d %d \n",m,rdex,idd,idk,ctg_used[rdex],cono,mono,ctg_part1[rdex],ctg_part2[rdex]);
                        if((ctg_length[rdex] < min_len)&&(ctg_used[rdex] == 0))
                        {
                      printf("hhh0: %d %d %d %d %d | %d %d || %d %d \n",m,rdex,idd,idk,ctg_used[rdex],cono,mono,ctg_part1[rdex],PO_matrix[idd][m]);
                          if((PO_matrix[idd][m]%2) == 0)
                            mono = 2;
                          else
                            mono = 1;
                          if(cono != mono)
                          {
                            min_len = ctg_length[rdex];
                            min_idd = m;
                            idk = PN_matrix[idd][m];
                            if((PO_matrix[idd][m] >= 2)&&(PO_matrix[idd][m] <= 3))
                              rc_idk = 0;
                            else
                              rc_idk = 1;
                          } 
                        }
                      printf("hhh1: %d %d %d %d %d %d %d || %d %d \n",m,rdex,idd,idk,ctg_used[rdex],ctg_used[idk],ptn_list[idd],ctg_part1[rdex],ctg_part2[rdex]);
                     }
                   }
                 }
                 printf("hhhx: %d %d %d %d %d %d | %d %d\n",k,idd,idk,ctg_part1[idd],ctg_part2[idd],rc_idk,ctg_used[idd],ctg_used[idk]);
                 if((ptn_list[k]==1)&&(ctg_used[k]==0)&&(ctg_print[k] == 1))
                 {
                   printf("supercontig2: tarseq_%d %d %d %d %d %d\n",k,idd,idk,rc_idk,ctg_rcdex1[k],ctg_patnum[k]);
                   if(rc_idk==0)
                   {
                     if(rcdex_matrix[k][idk] == 1)
                     {
                       fprintf(namef,"contigg1:%08d %d %d %d 1 2\n",0,n_scaff,k,ctg_length[k]);
                       printf("contigg1: %d %d %d 0 %d\n",n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                     }
                     else if(rcdex_matrix[k][idk] == 3)
                     {
                       fprintf(namef,"contigg1:%08d %d %d %d 0 3\n",0,n_scaff,k,ctg_length[k]);
                       printf("contigg1: %d %d %d 1 %d\n",n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                     }
                     else
                     {
                       fprintf(namef,"contigg1:%08d %d %d %d 0 %d\n",0,n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                       printf("contigg1: %d %d %d 0 %d\n",n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                     }
                     ctg_output[k] = 1;
                     ctg_outrc[k] = 0;
                     if(idk!=k)
                     {
                       if((ctg_output[idk] == 0)&&(idk >= 0))
                       {
                         fprintf(namef,"contigg2:%08d %d %d %d 0 0\n",s2_matrix[k][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 0 0\n",n_scaff,idk,ctg_length[idk]);
                         ctg_outrc[idk] = 0;
                         ctg_part2[k] = idk;
                         ctg_part1[idk] = k;
                       }
                       ctg_output[idk] = 1;
                     }
                   }
                   else if(rc_idk==1)
                   {
                     fprintf(namef,"contigg1:%08d %d %d %d 0 %d\n",0,n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                     printf("contigg1: %d %d %d 0 %d\n",n_scaff,k,ctg_length[k],rcdex_matrix[k][idk]);
                     ctg_output[k] = 1;
                     ctg_outrc[k] = 0;
                     if(idk!=k)
                     {
                       if((ctg_output[idk] == 0)&&(idk >= 0))
                       {
                         fprintf(namef,"contigg2:%08d %d %d %d 1 0\n",s2_matrix[k][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 1 0\n",n_scaff,idk,ctg_length[idk]);
                         ctg_outrc[idk] = 1;
                         ctg_part2[k] = idk;
                         ctg_part1[idk] = k;
                       }
                       ctg_output[idk] = 1;
                     }
                   }
                   else
                   {
                     fprintf(namef,"contigg1:%08d %d %d %d 0 4\n",0,n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 0 0\n",n_scaff,k,ctg_length[k]);
                     ctg_output[k] = 1;
                   }
                   tbase = tbase + ctg_length[k];
                   if(idk!=k)
                     tbase = tbase + ctg_length[idk];
                   ctg_used[k] = 0;
                   num_loops++;
                 }
                 else if((ctg_print[idd]==1)&&(idd!=idk))
                 {
                   int rc_idd = 0;
                   int rc_ide = 0;

//                   if(idk>=0)
                   tbase = tbase + ctg_length[idk];
                   if(rc_idk==0)
                     rc_idd = 0;
                   else if(rc_idk==1)
                     rc_idd = 1;
                   else if(rc_idk==2)
                     rc_idd = 0;
                   else if(rc_idk==3)
                     rc_idd = 1;
                   if(rc_idd!=rc_idi)
                     rc_ide = 1;
                   else
                     rc_ide = 0;
                   if(ctg_output[idk] == 0)
                   {
                     int outrc = 0;
                     if(ctg_outrc[idd] == 0)
                       outrc = rc_matrix[idd][idk];
                     else
                     {
                       if(rc_matrix[idd][idk] == 0)
                         outrc = 1;
                       else
                         outrc = 0;
                     }
                     ctg_outrc[idk] = outrc;
                     if(idk >= 0)
                     {
                       fprintf(namef,"contig-0:%08d %d %d %d %d 0\n",s2_matrix[idd][idk],n_scaff,idk,ctg_length[idk],outrc);
/*
                       fprintf(namef,"         PP: %d %d %d ",n_scaff,idk,ptn_list[idk]);
                       for(m=0;m<ptn_list[idk];m++)
                          fprintf(namef,"%d ",PN_matrix[idk][m]);
                       fprintf(namef,"|| ");
                       for(m=0;m<ptn_list[idk];m++)
                          fprintf(namef,"%d ",PO_matrix[idk][m]);
                       fprintf(namef,"\n");
            */
                       printf("contig-0: %d %d %d %d || %d %d %d %d %d\n",n_scaff,idk,ctg_length[idk],outrc,idd,idk,ctg_outrc[idd],ctg_outrc[idk],rc_matrix[idd][idk]);
                       ctg_part2[idd] = idk;
                       ctg_part1[idk] = idd;
                     }
                   }
                   ctg_output[idk] = 1;
                   rc_idi = rc_ide;       
                   num_loops++;
                 }
//             if(idd==i_getindex)
               printf("hhh3: %d %d %d\n",k,idd,idk);
                 ctg_used[k] = 1;
                 ctg_used[idd] = 1;
                 idd = idk;
               }
               else
                 stopflag=1;
               if(stopflag == 1)
                 break;
               printf("hhh-xxx: %d %d %d %d %d %d %d\n",k,idd,idk,stopflag,ctg_print[idd],ctg_links[idk],ctg_mask[idk]);
             }
             if(tbase == 0)
               tbase = ctg_length[idk];
             if(num_loops != 0)
               printf("bases: %d %d %d\n",idk,n_scaff,tbase);
             n_scaff++;
            }
          }
          for(k=0;k<max_ctg;k++)
          {
             int tbase = 0;
             if(ctg_output[k] == 0)
             {
               printf("supercontig3: tarseq_%d %d %d\n",k,k,ctg_patnum[k]);
               fprintf(namef,"contig-n:%08d %d %d %d 0 0\n",0,n_scaff,k,ctg_length[k]);
               printf("contig-n: %d %d %d 0 0\n",n_scaff,k,ctg_length[k]);
               tbase = ctg_length[k];
               printf("bases: %d %d %d\n",k,n_scaff,tbase);
               n_scaff++;
             }
          }
     fclose(namef);
}

/* ========================================================= */
void Reads_Distribution(int nSeq, int R_len, int *rd_locus, int *rd_index)
/* ========================================================= */
{
    int i,j,k,stopflag,num_steps,num_ave,*n_hit;
    int num_hits,BAR = 5000;
    int nstep = 20000;
    double rate,rate2;

    if((n_hit= (int *)calloc(nSeq,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: n_hit\n");
      exit(1);
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
         if((rd_locus[j]<(BAR+nstep))&&(rd_locus[i]>=BAR))
         {
           j++;
         }
         else
           stopflag=1;
       }
       if((j-i)>=3)
       {
         rate = (j-i)*100;
         rate = rate/nSeq;
//         if((i_getindex == 3)&&(j_getindex == 52))
//           printf("frequency:%d %d %f\n",j-i,BAR,rate);
         n_hit[num_hits] = j-i;
         BAR = BAR+nstep;
         num_hits++;
         num_steps++;
       }
       else if((j-i)<=2)
       {
         rate = 100;
         rate = rate/nSeq;
         BAR = rd_locus[i];
//         num_steps++;
       }
       i=j-1;
     }

     rate2 = 0.0;
     rate = R_len;
     rate = rate/nstep;
     if(num_steps > 0)
     {
       rate2 = num_steps;
       rate2 = rate2/rate;
     }
     else
       rate2 = 0.0;

     c_score = rate2;

     num_ave = nSeq;
     if(num_hits == 0)
       c_pairs = nSeq;
     else
     {
       num_ave = num_ave/num_hits;
       num_ave = num_ave*5;
       c_pairs = 0;
       for(i=0;i<num_hits;i++)
       {
          if(n_hit[i] <= num_ave)
            c_pairs = c_pairs+n_hit[i];
       }
     }
//     printf("Num_steps: %d %lf %lf %d\n",num_steps,rate,rate2,R_len);

}


/* ======================================= */
void Direction_Ratio(int nSeq,int *dlinks) 
/* ======================================= */
{
     double mf1,mf2,sigama;
     int dlink_index[5];
     int i,j,k;
     void ArraySort2_Int2(int n, int *arr, int *brr);

     for(i=0;i<4;i++)
        dlink_index[i] = i;
     ArraySort2_Int2(4,dlinks,dlink_index);

     mf1 = dlinks[1];
     mf1 = mf1/dlinks[0];
     if(mf1 > 0.9)
       dlinks[4] = 0;
     else
       dlinks[4] = 1;
}


/* ========================================================= */
void Normal_Distribution(int nSeq, int R_len, int *rd_locus, int *rd_index)
/* ========================================================= */
{
     int i,j,k,stopflag,num_steps,num_hits,num_reads;
     int hit_max,hit_loc,hit_buk,hit_siga;
     int *n_hit,*s_len,BAR = 5000;
     int nstep = 20000;
     double rate,rate2,sigama;

     if((n_hit= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: n_hit\n");
       exit(1);
     }
     if((s_len= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: s_len\n");
       exit(1);
     }

     num_reads = 0;
     num_hits = 0;
     num_steps = 0; 
     for(i=0;i<nSeq;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if((rd_locus[j]<(BAR+nstep))&&(rd_locus[i]>=BAR))
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>1)
        {
          rate = (j-i)*100;
          rate = rate/nSeq;
          BAR = BAR+nstep;
          n_hit[num_hits] = j-i;
          s_len[num_hits] = rd_locus[i];
//      printf("Array: %d %d %d %d %lf\n",nSeq,num_hits,j-i,s_len[num_hits],rate);
          num_hits++;
          num_steps++;
        }
        else if((j-i)==1)
        {
          rate = 100;
          rate = rate/nSeq;
          BAR = rd_locus[i];
        }
        i=j-1;
     }

     hit_max = 0;
     hit_loc = 0;
     hit_buk = 0;

     nSeq = num_hits;
     for(i=0;i<nSeq;i++)
     {
        if(n_hit[i] > hit_max)
        {
          hit_max = n_hit[i];
          hit_loc = s_len[i];
          hit_buk = i;
        }
        num_reads = num_reads + n_hit[i];
     }

    hit_siga = 0;
    hit_max = 0;
    rate = 0.0;
    for(i=hit_buk;i<nSeq;i++)
    {
       hit_max = hit_max + n_hit[i];
       sigama = hit_max;
       sigama = sigama/num_reads;
       if(sigama >= 0.341)
       {
         hit_siga = s_len[i]-hit_loc;
         rate = s_len[i]-hit_loc;
         rate =  rate/s_len[nSeq-1];
//         printf("AA: %d %d %lf\n",i,s_len[i],rate);
         i = nSeq;
       } 
    }

/*
    if(rate == 0.0)
      printf("F: %d %d %lf %d\n",i,s_len[nSeq-1],rate,nSeq);
    else
      printf("F: %d %d %lf\n",i,s_len[i],rate);   */

    NL_score = rate;

    hit_max = 0;
    rate = 0.0;
    j = 0;
    for(i=hit_buk;i>=0;i--)
    {
       hit_max = hit_max + n_hit[i];
       sigama = hit_max;
       sigama = sigama/num_reads;
       if(sigama >= 0.341)
       {

         hit_siga = s_len[i]-hit_loc;
         rate = hit_loc-s_len[i];
         rate =  rate/s_len[nSeq-1];
//         printf("BB: %d %d %lf\n",i,s_len[i],rate);
         j = i;
         i = -1;
       } 
    }

/*
    if(rate == 0.0)
      printf("B: %d %d %lf\n",0,s_len[0],rate);
    else
      printf("B: %d %d %lf\n",j,s_len[j],rate);  */

    NR_score = rate;

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


