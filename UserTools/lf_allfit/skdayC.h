/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_INC_SKDAY_H
#define FH2H_INC_SKDAY_H


#ifdef __cplusplus
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif


/*------ fortran header (without commons and data statements) ----------*/

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*	skday.h*/
/**/
/*	definition of common variable of run*/
/**/
/*	(data and author)	97/3/17  H.Ishino*/
/*                               00/4/30  msm  (modified to prevent R10K)*/
/*                               05/10/15 Y.Takeuchi  max_run = 30000 -> 50000*/
/*                               08/05/22 Y.Takeuchi  max_run = 50000 -> 70000*/
/*                               12/07/17 Y.Takeuchi  max_run = 70000 ->999999*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
#ifndef SKDAY_H
#define SKDAY_H
#define MAX_RUN (999999)
/*year of the run*/
/*month of the run*/
/*day of the run*/
/*hour of the run*/
/*minites of the run*/
/*elpased day from 96/1/1*/
/*run time of the run*/

/*common skday_data was here*/
#endif 


/*------ common blocks -------------------------------------------------*/

extern struct skday_data_common {
  int    ryear[(MAX_RUN)-(1000)+1];
  int    rmon[(MAX_RUN)-(1000)+1];
  int    rday[(MAX_RUN)-(1000)+1];
  int    rhour[(MAX_RUN)-(1000)+1];
  int    rmin[(MAX_RUN)-(1000)+1];
  float  rlive[(MAX_RUN)-(1000)+1];
  int    relapse[(MAX_RUN)-(1000)+1];
} skday_data_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct skday_data_common *skday_data;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct skday_data_common *skday_data = &skday_data_;
#endif


/*------ data statements -----------------------------------------------*/


#ifndef NO_STATIC_DATA


#endif  /* #ifndef NO_STATIC_DATA */


/*------ end of fortran header -----------------------------------------*/


#ifdef __cplusplus
}
#endif


#endif  /* #ifndef FH2H_INC_SKDAY_H */
