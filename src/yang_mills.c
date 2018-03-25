#ifndef YM_VEC_C
#define YM_VEC_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#ifndef ONE_FILE_MODE
  #include"../include/function_pointers.h"
  #include"../include/gauge_conf.h"
  #include"../include/geometry.h"
  #include"../include/gparam.h"
  #include"../include/random.h"
#else
  #include"../include/endianness.h"
  #include"../include/function_pointers.h"
  #include"../include/gauge_conf.h"
  #include"../include/geometry.h"
  #include"../include/gparam.h"
  #include"../include/macro.h"
  #include"../include/mymalloc.h"
  #include"../include/random.h"
  #include"../include/su2.h"
  #include"../include/su2_upd.h"
  #include"../include/sun_aux.h"
  #include"../include/sun.h"
  #include"../include/sun_upd.h"

  #include"../lib/endianness.c"
  #include"../lib/function_pointers.c"
  #include"../lib/gauge_conf_def.c"
  #include"../lib/gauge_conf_meas.c"
  #include"../lib/gauge_conf_transfer.c"
  #include"../lib/gauge_conf_upd.c"
  #include"../lib/geometry.c"
  #include"../lib/gparam.c"
  #include"../lib/mymalloc.c"
  #include"../lib/random.c"
  #include"../lib/su2.c"
  #include"../lib/su2_upd.c"
  #include"../lib/sun_aux.c"
  #include"../lib/sun.c"
  #include"../lib/sun_upd.c"
#endif


void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
    int count;
    FILE *datafilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize function_pointers
    init_function_pointers();

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update_w(&GC, &geo, &param);

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures_localobs(&GC, &geo, &param, datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           save_on_file(&GC, &param);

           // backup copy
           save_on_file_back(&GC, &param);
           }
         }

       // save configuration for offline analysis
       if(param.d_saveconf_analysis_every!=0)
         {
         if(count % param.d_saveconf_analysis_every == 0 )
           {
           strcpy(name, param.d_conf_file);
           sprintf(aux, "%ld", GC.update_index);
           strcat(name, aux);
           save_on_file_with_name(&GC, &param, name);
           }
         }
       }
    time(&time2);
    // montecarlo end

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      save_on_file(&GC, &param);
      }

    // print simulation details
    print_parameters(&param, time1, time2);

    // free gauge configuration
    end_gauge_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);

    exit(EXIT_SUCCESS);
    }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("Program %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef ONE_FILE_MODE
        printf("\n\tcompiled in the single file mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

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
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
