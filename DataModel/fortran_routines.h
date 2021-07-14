// (c-versions of) headers that define fortran common blocks
// these are offficial in $SKOFL_ROOT/inc/
#include "skheadC.h"
#include "skparmC.h"
#include "sktqC.h"
#include "skbadcC.h"
#include "geopmtC.h"
#include "skruninfC.h"  // commented out in original lowfit_sk4

// these are from $SKOFL_ROOT/inc/lowe and have no official C version
// the following two were done via automatic conversion done with fh2h
#include "skdayC.h"
#include "skwtC.h"
// this last one was manually edited to reflect 'EQUIVALENCE' use (use with caution)
#include "skroot_loweC.h"

// header for skroot_* functions. These are actually C functions.
#include "fortran_interface.h"

// we need to declare all external fortran routines used
// the following were fine without any additional libraries being added to the makefile
extern "C" void skroot_init_(int*);
extern "C" void kzinit_();
extern "C" void skoptn_(char*, int);
extern "C" void skbadopt_(int*);
extern "C" void geoset_();
extern "C" void delete_outside_hits_();
extern "C" void skcrawread_(int*, int*);
extern "C" void skcread_(int*, int*);
extern "C" void skroot_set_tree_(int*);
extern "C" void skroot_get_entry_(int*);

// from $ATMPD_ROOT/src/programs/TreeBuilder/examples/fort_fopen.F
extern "C" void fort_fopen_(int*, const char*, char*, int* ,int);

// read runinf (y/m/d, start time, end time, etc)
extern "C" void runinfsk_();

// the following are provided by libwtlib_5.1.a
extern "C" void skrunday_();
extern "C" void skwt_();
extern "C" void skwt_gain_corr_();
extern "C" void lfwater_(int*, float*);
// skday_data_, common block

// the following are provided by libbonsai_3.3.a
extern "C" void cfbsinit_(int*, float*);
extern "C" void cfbsexit_();

// the following are provided by libsklowe_7.0.a
extern "C" void lfclear_all_();
extern "C" void lfallfit_sk4_final_qe43_(float*, int*, int*, int*, int*);
extern "C" void lfallfit_sk4_data_(float*, int*, int*);
extern "C" void lfallfit_sk4_gain_corr_(float*, int*, int*, int*, int*);
extern "C" void lfallfit_sk4_mc_(float*, int*, int*);
// skroot_lowe_ common block

// TODO where from
extern "C" void skbadch_(int*, int*, int*);

// after that there were many undefined references, e.g. `sortzv_`, `hf1_`, `hf2_`...
// after some trial and error these are resolved, but i lost track of which provided what.
// cernlibs in particular resolved a lot of repeated undefined issues, they may be the main culprit.
