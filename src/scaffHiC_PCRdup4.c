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
static int *hit_rddex,*hit_score,*hit_rcdex,*hit_locus1,*superlength,*hit_matlocu1,*hit_matlocu2,*hit_matindex;
static int *ctg_length,*hit_index,*hit_masks,*hit_linkdex1,*hit_linkdex2;;

/* SSAS default parameters   */
static int IMOD=0;
static int nContig=0;
static int file_flag = 1;
static int min_len = 100000;
static int max_ctg = 1000000;
static int n_depth = 60;
static float m_score = 400.0;
static int i_getindex = 2;

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
      printf("Usage: %s <input_readplace_file> <output_readplace_file>\n",argv[0]);

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
    if((hit_masks = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_masks\n");
      exit(1);
    }

    nSeq=nseq;
    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    S_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
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

    n_reads=i;
    Matrix_Process(argv,args,n_reads);

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
//        printf("%d %d %s",nseq,hit_score[nseq],line);
      if(hit_score[nseq] == 0)
        printf("%s",line);
      nseq++;
    }
    fclose(namef); 

    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Matrix_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,n_scaff;
     FILE *namef,*namef2;
     int *ctg_list;
     long num_cells,n_Bases,ii,n_blocks;
     int num_hits,num_hit1,num_hit2,rcdex,rsize,rsize2,size_row;
     int stopflag,offset,*ray,*dex1,*dex2;
     void ArraySort_Mix(int n, long *arr, int *brr);
     void ArraySort_float2(int n, float *arr, int *brr);
     char **DBname,*st,*ed,line[2000];
     void PCRdup_Process(char **argv,int args,int nSeq,int offset1,int offset2);
     int **p_matrix,**s_matrix,**s2_matrix,**o_matrix,**r_matrix,**rc0_matrix,**rc1_matrix,**rc2_matrix,**rc3_matrix,**rc_matrix;
     float rate,*Dis_index,*Dis_ratia1,*Dis_ratia2;
     int **imatrix(long nrl,long nrh,long ncl,long nch);
     float **fmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     void ArraySort_Int2(int n, int *arr, int *brr);
     int *ctg_score1,*ctg_score2,*ctg_mapp1,*ctg_mapp2,*ctg_join1,*ctg_join2,*ctg_idex1,*ctg_idex2,*ctg_mask,*ctg_rcdex1;
     int *p_index,*p_rcdex,*p_masks,*p_score,*p_lists;
     int *ctg_rcoo1,*ctg_rcoo2,*ctg_part1,*ctg_part2,*ctg_patnum;
     int *ctg_output,*ctg_hitnum,*ctg_used,*ctg_links,*ctg_oodex1,*ctg_oodex2,*ctg_rcindex,*ctg_mpindex,*ctg_outrc;
     int *link_locus,*link_locu2,*link_locu3,*head_locus;
     int n_length = min_len;

     rsize = max_ctg+10; 
     n_blocks = rsize*rsize; 
     printf("Memory: %ld %d %d\n",n_blocks,max_ctg,rsize);
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

     if((link_locus = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - link_locus\n");
       exit(1);
     }
     if((head_locus = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - head_locus\n");
       exit(1);
     }
     if((link_locu2 = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - link_locu2\n");
       exit(1);
     }
     if((link_locu3 = (int *)calloc(n_blocks,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - link_locu3\n");
       exit(1);
     }

     printf("contigs: %d %d %d\n",rsize,max_ctg,nSeq);
     size_row = n_depth;
     p_matrix=imatrix(0,rsize,0,rsize);
     s_matrix=imatrix(0,rsize,0,rsize);
     s2_matrix=imatrix(0,rsize,0,rsize);
     o_matrix=imatrix(0,rsize,0,rsize);
     r_matrix=imatrix(0,rsize,0,rsize);
     rc0_matrix=imatrix(0,rsize,0,rsize);
     rc1_matrix=imatrix(0,rsize,0,rsize);
     rc2_matrix=imatrix(0,rsize,0,rsize);
     rc3_matrix=imatrix(0,rsize,0,rsize);
     rc_matrix=imatrix(0,rsize,0,rsize);
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

     printf("contigs2: %d %d %d\n",rsize,max_ctg,nSeq);
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
     for(ii=1;ii<n_blocks;ii++)
     {
        head_locus[ii] = head_locus[ii-1] + link_locu3[ii-1];
     }

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
     if((hit_linkdex1 = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_linkdex1\n");
       exit(1);
     }
     if((hit_linkdex2 = (int *)calloc(num_cells,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_linkdex2\n");
       exit(1);
     }

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
          long idh1,idh2;

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
          hit_matlocu2[head_locus[idh2]+link_locu2[idh2]] = loci2;
          strcpy(S_Name[head_locus[idh1]+link_locus[idh1]],R_Name[i]);
          hit_linkdex1[head_locus[idh1]+link_locus[idh1]] = i;
          hit_linkdex2[head_locus[idh2]+link_locu2[idh2]] = i;
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
     for(i=0;i<max_ctg;i++)
     {
        for(j=0;j<max_ctg;j++)
           hit_rddex[j] = j;
        Dis_index[i] = 0.0;
        ArraySort2_Int2(max_ctg,r_matrix[i],hit_rddex);
 
        memset(ctg_rcindex,0,4*max_ctg); 

        printf("scaffold: %d %d %d\n",i,max_ctg,ctg_length[i]);

        printf("matrix: %d %d\n",i,ctg_length[i]);
        for(j=0;j<size_row;j++)
        {
           int idi,idj,offset1,offset2;
           int idd = hit_rddex[j];
           int n_pairs = r_matrix[i][j];
           int n_pairs_half = n_pairs/2;
           float rr1,rr2,rat1,rat2,sq1,sq2;
           int halflen_1,halflen_2,half_len1,half_len2;
           int rcindex1,rcindex2,rcindex;
   
           rcindex1 = 1;
           rcindex2 = 1;

           idi = i;
           idj = hit_rddex[j]; 
           
           offset1 = head_locus[idi*max_ctg+idj];
           ray = hit_matlocu1;
           dex1 = hit_linkdex1; 
           ArraySort_Int2(n_pairs,ray+offset1,dex1+offset1);

           offset2 = head_locus[idj*max_ctg+idi];
           ray = hit_matlocu2;
           dex2 = hit_linkdex2; 
           ArraySort_Int2(n_pairs,ray+offset2,dex2+offset2);


           if((i == 82)&&(idd == 86))
           {
             PCRdup_Process(argv,args,n_pairs,offset1,offset2);
           }
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
           rr1 = n_Bases;
           rr1 = rr1/nSeq;
           rr1 = rr1*n_pairs*1000.0; 
           rr1 = rr1/half_len1;
           rr1 = rr1*(rat1 - 1.0);
//           rr1 = rr1/half_len1;

           rr2 = n_Bases;
           rr2 = rr2/nSeq;
           rr2 = rr2*n_pairs*1000.0; 
           rr2 = rr2/half_len2;
           rr2 = rr2*(rat2 - 1.0);
//           rr2 = rr2/half_len2;

           if(i==1)
           {
//                printf("mmm: %d %ld %d %d %d %d %d %f %f %f %f\n",n_pairs,n_Bases,nSeq,halflen_2,ctg_length[idj]/2,idi,idj,rr1,rr2,sq1,sq2);
           }

           Dis_ratia1[j] = rat1;
           Dis_ratia2[j] = rat2;

           Dis_index[j] = rr1*rr2; 
//           printf("%6d | %6d %d %d | %d %d %f %f | %f %f %f %4d %4d %4d %4d\n",ctg_length[idd],r_matrix[i][j],hit_rddex[j],rcindex,halflen_1,halflen_2,rat1,rat2,rr1,rr2,rr1*rr2,rc0_matrix[i][idd],rc1_matrix[i][idd],rc2_matrix[i][idd],rc3_matrix[i][idd]);
        }

     }

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
     if((namef2 = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
     i=0;
     while(!feof(namef))
     {
       if(fgets(line,2000,namef) == NULL)
       {
//        printf("fgets command error:\n);
       }
       if(feof(namef)) break;
       if(hit_masks[i] == 0)
         fprintf(namef2,"%s",line);
       i++;
     }
     fclose(namef);
     fclose(namef2);
     printf("Masked reads %d\n",i);
}

/*   subroutine to sort out read pairs    */
/* =============================== */
void PCRdup_Process(char **argv,int args,int nSeq,int offset1,int offset2)
/* =============================== */
{
     int i,j,k,m,n,n_scaff,n_blocks;
     int BAR,stopflag,nstep,base,num_steps;
     float rate;
     int *heap_list1,*heap_list2,*heap_steps,num_heap1,num_heap2;

     if((heap_list1 = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }
     if((heap_list2 = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }
     if((heap_steps = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }

     num_heap1 = 0;
     num_heap2 = 0;
     BAR = 5000;
     nstep = 20000;
     num_steps = 0;

//     for(k=0;k<nSeq;k++)
//        printf("kkk: %d %d %d %s\n",k,hit_matlocu1[offset1+k],hit_matlocu1[offset2+k],R_Name[hit_linkdex2[offset2+k]]);                
     for(k=0;k<nSeq;k++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=k+1;
        base = hit_matlocu1[offset1+k];; 
        while((j<nSeq)&&(stopflag==0))
        {
          if((hit_matlocu1[offset1+j]<(BAR+nstep))&&(hit_matlocu1[offset1+k]>=BAR))
          {
            base = base + hit_matlocu1[offset1+j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-k)>1)
        {
          rate = (j-k)*100;
          rate = rate/nSeq;
          printf("frequency1:%d %d %f\n",j-k,BAR,rate);
          BAR = BAR+nstep;
          heap_list1[num_heap1] = j-k;
          heap_steps[num_heap1] = k; 
          num_heap1++;
          num_steps++;
        }
        else if((j-k)==1)
        {
          rate = 100;
          rate = rate/nSeq;
//          printf("frequency11: %d %f\n",BAR,rate);
          BAR = hit_matlocu1[offset1+k];
          num_steps++;
        }
        k = j-1;
     }

     num_hits = 0;
     for(i=0;i<num_heap1;i++)
         num_hits = num_hits + heap_list1[i];
    
     num_ave = num_hits;
     if(num_heap1 > 0)
       num_ave = num_hits/num_heap1;
     else
       num_ave = num_hits;

     for(i=0;i<num_heap1;i++)
     {
        int rat1,rat2;

        if(i==0)
        {
          if(heap_list1[i+1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i+1];
          rat2 = rat1;
        }
        else if(i==(num_heap1-1))
        {
          if(heap_list1[i-1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i-1];
          rat2 = rat1;
        }
        else
        {
          rat1 = heap_list1[i]/heap_list1[i-1];
          rat2 = heap_list1[i]/heap_list1[i+1];
        }
        if((rat1 >= 3)&&(rat2 >=3))
        {
          int num_dups = heap_list1[i];
          int idk = heap_steps[i];
          for(j=0;j<num_dups;j++)
          {
             int idd = hit_linkdex1[offset1+idk+j];
             if((idd%2) == 0)
             {
               hit_masks[idd] = 1;
               hit_masks[idd+1] = 1;
             }
             else
             {
               hit_masks[idd] = 1;
               hit_masks[idd-1] = 1;
             }
          }
        }
     }
 //    printf("hit1: %d %d\n",num_heap1,nSeq);


     num_heap1 = 0;
     num_heap2 = 0;
     BAR = 5000;
     nstep = 20000;
     num_steps = 0;
     base = 0;
     memset(heap_list1,0,4*nSeq);
     memset(heap_steps,0,4*nSeq);

//     for(k=0;k<nSeq;k++)
//        printf("kkk: %d %d %d %s\n",k,hit_matlocu1[offset1+k],hit_matlocu2[offset2+k],R_Name[hit_linkdex2[offset2+k]]);

     for(k=0;k<nSeq;k++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=k+1;
        base = hit_matlocu2[offset2+k];; 
        while((j<nSeq)&&(stopflag==0))
        {
          if((hit_matlocu2[offset2+j]<(BAR+nstep))&&(hit_matlocu2[offset2+k]>=BAR))
          {
            base = base + hit_matlocu2[offset2+j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-k)>1)
        {
          rate = (j-k)*100;
          rate = rate/nSeq;
          printf("frequency2:%d %d %f\n",j-k,BAR,rate);
          BAR = BAR+nstep;
          heap_list1[num_heap1] = j-k;
          heap_steps[num_heap1] = k; 
          num_heap1++;
          num_steps++;
        }
        else if((j-k)==1)
        {
          rate = 100;
          rate = rate/nSeq;
//          printf("frequency22: %d %f\n",BAR,rate);
          BAR = hit_matlocu2[offset2+k];
          num_steps++;
        }
        k = j-1;
     }

     for(i=0;i<num_heap1;i++)
     {
        int rat1,rat2;

        if(i==0)
        {
          if(heap_list1[i+1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i+1];
          rat2 = rat1;
        }
        else if(i==(num_heap1-1))
        {
          if(heap_list1[i-1] == 0)
            rat1 = 10; 
          else
            rat1 = 3*heap_list1[i]/heap_list1[i-1];
          rat2 = rat1;
        }
        else
        {
          rat1 = heap_list1[i]/heap_list1[i-1];
          rat2 = heap_list1[i]/heap_list1[i+1];
        }
        if((rat1 >= 3)&&(rat2 >=3))
        {
          int num_dups = heap_list1[i];
          int idk = heap_steps[i];
          for(j=0;j<num_dups;j++)
          {
             int idd = hit_linkdex2[offset2+idk+j];
             if((idd%2) == 0)
             {
               hit_masks[idd] = 1;
               hit_masks[idd+1] = 1;
             }
             else
             {
               hit_masks[idd] = 1;
               hit_masks[idd-1] = 1;
             }
          }
        }
     }
//     printf("hit2: %d %d\n",num_heap1,nSeq);
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


