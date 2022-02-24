#ifndef MEASURE_POLY_C
#define MEASURE_POLY_C

#include"../include/macro.h"
#include<openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/function_pointers.h"
#include"../include/gparam.h"
#include"../include/gauge_conf.h"



int main(int argc, char **argv)
   {
   char infile[STD_STRING_LENGTH];
   GParam param;
   Gauge_Conf GC;
   Geometry geo;
   FILE*fp, *datafilep, *monofilep, *cluster_monofilep, *coords;
   int dim, err, i;

   if(argc != 2)
     {
     printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
     printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
     printf("Usage: %s conf_file\n\n", argv[0]);

     printf("Compilation details:\n");
     printf("\tGGROUP (gauge group): %s\n", QUOTEME(GGROUP));
     printf("\tN_c (number of colors): %d\n", NCOLOR);
     printf("\n");
     printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
     printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

     return EXIT_SUCCESS;
     }
   else
     {
     if(strlen(argv[1]) >= STD_STRING_LENGTH)
       {
       fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
       }
     else
       {
       strcpy(infile, argv[1]);
       }
     }

   //Open the configuration file and save the dimensions in param
   fp=fopen(infile, "r");
   if(fp==NULL)
     {
     fprintf(stderr, "Error in opening the file %s (%s, %d)\n", infile, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   else
     {
     err=fscanf(fp, "%d", &dim);
     if(err!=1)
       {
       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", infile, __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     for(i=0; i<dim; i++)
        {
        err=fscanf(fp, "%d", &(param.d_size[i]));
        if(err!=1)
          {
          fprintf(stderr, "Error in reading the file %s (%s, %d)\n", infile, __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

   fclose(fp);

   // things to be initialized in order to read the conf
   param.d_start=2;
   strcpy(param.d_conf_file, infile);
   init_derived_constants(&param);
   init_indexing_lexeo();
   init_geometry(&geo, &param);
   init_gauge_conf(&GC, &param);

   // open the ouputfile
   datafilep = fopen("plaq_pol.dat", "a");
   monofilep = fopen("mon.dat", "a");
   cluster_monofilep = fopen("mon_clusters.dat", "a");
   coords = fopen("coords.dat", "a");

   // compute the Polyakov observables
   double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1], mod_loc, aux;

   plaquette(&GC, &geo, &param, &plaqs, &plaqt);
   polyakov_for_tracedef(&GC, &geo, &param, polyre, polyim);
   polyakov_loc_fluct(&GC, &geo, &param, polyre, polyim, &mod_loc);
   fprintf(datafilep, "%ld %.12g %.12g ", GC.update_index, plaqs, plaqt);

   for(i=0; i<(int)floor(NCOLOR/2); i++)
      {
      fprintf(datafilep, "%.12g %.12g %.12g %.12g ", polyre[i], polyim[i], NCOLOR*sqrt(pow(polyre[i],2)+pow(polyim[i],2)),
                                                           mod_loc - NCOLOR*sqrt(pow(polyre[i],2) + pow(polyim[i],2)));
      }
   aux = NCOLOR*sqrt(pow(polyre[0],2) + pow(polyim[0],2));
   fprintf(datafilep, "\n");   

   fclose(datafilep);

   Gauge_Conf helperconf;
   int subg, subgnum;

   init_gauge_conf_from_gauge_conf(&helperconf, &GC, &param);
   alloc_diag_proj_stuff(&helperconf, &param);

   // MAG gauge fixing
   max_abelian_gauge_fix(&helperconf, &geo, &param);

   //diagonal projection
   diag_projection(&helperconf, &param);

   //loop on all the U(1) subgroups
   if(NCOLOR>1)
     {
     subgnum=NCOLOR-1;
     }
   else
     {
     subgnum=1;
     }
   for(subg=0; subg<subgnum; subg++)
      {
    
      // extract the abelian component subg and save it to GC->u1_subg

      U1_extract(&helperconf, &param, subg);
      // compute monopole observables
      monopoles_obs(&helperconf, &geo, &param, subg, aux, monofilep, coords);

      // compute the cluster of monopoles currents.
      monopoles_clusters(&helperconf, &geo, &param, subg, cluster_monofilep);

      }
      free_diag_proj_stuff(&helperconf, &param);
      free_gauge_conf(&helperconf, &param);

      fflush(monofilep);
      fclose(monofilep);
      
      fflush(coords);
      fclose(coords);
            
      free_gauge_conf(&GC, &param);
      free_geometry(&geo, &param);


      return 0;
      }
#endif

