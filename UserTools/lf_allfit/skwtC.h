/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_INC_SKWT_H
#define FH2H_INC_SKWT_H


#ifdef __cplusplus
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif


/*------ fortran header (without commons and data statements) ----------*/

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*	skwt.h*/
/**/
/*	defintion of variables for information of water transparancy and date*/
/**/
/*	(date and author)	97/3/17    H.Ishino*/
/*                               00/4/30    msm (modified to allow wt*/
/*                                               calculations past 1400 days)*/
/*                               05/10/15   Y.Takeuchi */
/*                                            wtday_num 400 -> 2000*/
/*                               08/05/22   Y.Takeuchi */
/*                                            wtday_num 2000 -> 8000*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

#ifndef SKWT_H
#define SKWT_H

#define WTDAY_NUM (8000)
/* number of day*/
/* max water transparancy*/
/* miniman water transparancy*/
/* day */
/* water transparancy(cm)*/

/*common wttable was here*/
/*common wttablenum was here*/
/*common wttablemm was here*/

#endif 


/*------ common blocks -------------------------------------------------*/

extern struct wttable_common {
  int    wtday[WTDAY_NUM];
  float  wtvalue[WTDAY_NUM];
} wttable_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct wttable_common *wttable;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct wttable_common *wttable = &wttable_;
#endif

extern struct wttablenum_common {
  int    wtday_exit_num;
} wttablenum_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct wttablenum_common *wttablenum;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct wttablenum_common *wttablenum = &wttablenum_;
#endif

extern struct wttablemm_common {
  float  wt_value_max;
  float  wt_value_min;
} wttablemm_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct wttablemm_common *wttablemm;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct wttablemm_common *wttablemm = &wttablemm_;
#endif


/*------ data statements -----------------------------------------------*/


#ifndef NO_STATIC_DATA


#endif  /* #ifndef NO_STATIC_DATA */


/*------ end of fortran header -----------------------------------------*/


#ifdef __cplusplus
}
#endif


#endif  /* #ifndef FH2H_INC_SKWT_H */
