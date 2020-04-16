#ifndef ROTATE_POLY_C
#define ROTATE_POLY_C

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
   FILE*fp, *datafilep;
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
   datafilep = fopen("plaq_pol.dat", "w");

   // compute the Polyakov observables
   double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1];

   plaquette(&GC, &geo, &param, &plaqs, &plaqt);
   polyakov_for_tracedef(&GC, &geo, &param, polyre, polyim);
   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

   for(i=0; i<(int)floor(NCOLOR/2); i++)
      {
      fprintf(datafilep, "%.12g %.12g %.12g ", polyre[i], polyim[i], NCOLOR*sqrt(pow(polyre[i],2)+pow(polyim[i],2)));
      }

   // now the rotation

   double angle;
   double tolerance = PI/3.0;
   double complex  z_pos, z_neg;

   z_pos = -0.5 +  sqrt(3)/2.0*I; // 2pi/3
   z_neg = -0.5 -  sqrt(3)/2.0*I; // 4pi/3

   angle = atan2(polyim[0], polyre[0]);

   //case 2pi/3
   if((angle< (2.0*PI/3)+tolerance) && (angle >(2.0*PI/3)-tolerance))
   {
   for(long r=0; r<param.d_space_vol;r++)
      {
      r = sisp_and_t_to_si(&geo, r, 0);
      for(int mu=0;mu<STDIM;mu++)
         {
         times_equal_complex(&(GC.lattice[r][mu]), z_pos);
         }
      }
   }
   // 4pi/3 case
   else if ((angle<(-2.0*PI/3)-tolerance) && (angle >-(2*PI/3)+tolerance))
   {
   for(long r=0; r<param.d_space_vol;r++)
      {
      r = sisp_and_t_to_si(&geo, r, 0);
      for(int mu=0;mu<STDIM;mu++)
         {
         times_equal_complex(&(GC.lattice[r][mu]), z_neg);
         }
      }
   }

   char aux[STD_STRING_LENGTH]="rot_";
   strcat(aux, param.d_conf_file);
   write_conf_on_file_with_name(&GC, &param, aux);

   plaquette(&GC, &geo, &param, &plaqs, &plaqt);
   polyakov_for_tracedef(&GC, &geo, &param, polyre, polyim);
   fprintf(datafilep, "%ld %.12g %.12g", GC.update_index, plaqs, plaqt);

   for(i=0; i<(int)floor(NCOLOR/2); i++)
      {
      fprintf(datafilep, "%.12g %.12g %.12g ", polyre[i], polyim[i], NCOLOR*sqrt(pow(polyre[i],2)+pow(polyim[i],2)));
      }


   fclose(datafilep);
   free_gauge_conf(&GC, &param);
   free_geometry(&geo, &param);


   return 0;
   }

#endif
