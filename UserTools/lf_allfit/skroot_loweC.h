/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_INC_SKROOT_LOWE_H
#define FH2H_INC_SKROOT_LOWE_H


#ifdef __cplusplus
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif


/*------ fortran header (without commons and data statements) ----------*/

/**/
/*  skroot_lowe.h     03-JUL-2008    Y.Takeuchi*/
/**/
/*  definition of variables for skroot_lowe files*/
/**/
/*  09-MAY-2008 spaevnum,spadll,spadl->spadlt, muboy_* added  by y.t*/
/*  22-MAY-2008 added SK-I parmeters, cleffwal, bseffwal*/
/*  08-MAR-2014 added Ariadne parameters*/
/*  07-OCT-2014 added T2K parameters by Koshio*/
/*  09-JUN-2016 added software trigger parameters by takeuchi*/
/**/
/* bonsai fit*/
/* clusfit*/
/* clusfit*/
/* time to the previous LE event (in raw data)*/
/* NS ratio*/
/* solar direction*/
/* event number of the parent muon*/
/* spallation log likelihood*/
/* longitudinal distance*/
/* traversal distance (usual delta-L)*/
/* MC vertex*/
/* MC direction (1st and 2nd particles)*/
/* MC absolute momentum (1st & 2nd)*/
/* MC energy (1st & 2nd)*/
/* MC dark rate for generation*/

      

/*common skroot_lowe was here*/

/**** muon info.*/
/* muboy status*/
/* number of tracks*/
/* up to 10 tracks*/
/* common direction*/
/* goodness*/
/* track length*/
/* dE/dX histogram*/
/* bff entpoint*/
/* bff direction*/
/* bff goodness*/
/* number of additional data in muinfo()*/
/* additional muon data (integer)*/
/* additional muon data (real)*/
/*EQUIVALENCE (muinfo(1),rmuinfo(1))*/
      

/*common skroot_mu was here*/

/**** SLE info.*/
/* bonsai fit*/
/* clusfit*/

/*common skroot_sle was here*/
    
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/**** muinfo ***************************************************************/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

/**** rmuinfo=1: original qismsk before fix_maxqisk()*/

/**** muinfo=2: number of sub-trigger */

/**** rmuinfo=3: original qimxsk before fix_maxqisk()*/

/**** rmuinfo=4: qimxsk after fix_maxqisk()*/
/*EQUIVALENCE (muinfo(4), muqimxsk)*/
      

/**** set muinfo=5: parent muon event number*/

/**** set muinfo=6: number of subtrigger in AFT*/

          
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/**** muinfo ***************************************************************/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/


/*------ common blocks -------------------------------------------------*/

extern struct skroot_lowe_common {
  double ltimediff;
  float  bsvertex[4];
  float  bsdir[3];
  float  bsresult[6];
  float  bsgood[3];
  float  bsdirks;
  float  bseffhit[12];
  float  bsenergy;
  int    bsn50;
  float  bscossun;
  float  clvertex[4];
  float  clresult[4];
  float  cldir[3];
  float  clgoodness;
  float  cldirks;
  float  cleffhit[12];
  float  clenergy;
  int    cln50;
  float  clcossun;
  int    latmnum;
  int    latmh;
  int    lmx24;
  float  lnsratio;
  float  lsdir[3];
  int    spaevnum;
  float  spaloglike;
  float  sparesq;
  float  spadt;
  float  spadll;
  float  spadlt;
  int    spamuyn;
  float  spamugdn;
  float  posmc[3];
  float  dirmc[2][3];
  float  pabsmc[2];
  float  energymc[2];
  float  darkmc;
  int    islekeep;
  float  bspatlik;
  float  clpatlik;
  float  lwatert;
  int    lninfo;
  union {
    int    linfo[255];
    float  rlinfo[255];   // /*EQUIVALENCE (linfo(1),rlinfo(1))*/
    /**** linfo ***************************************************************/
    struct {
        float nflf;                /**** linfo=1: flasher */
        float lnahit;              /**** linfo=2: lnahit*/
        float lnqisk;              /**** linfo=3: lnqisk*/
        float dist;                /**** linfo=4: dist*/
        float cleffwal;            /**** linfo=5: cleffwal*/
        float bseffwal;            /**** linfo=6: bseffwal*/
        float cleffh;              /**** linfo=7: cleffh    ! clusfit effective hit at lwatert*/
        float bseffh;              /**** linfo=8: bseffh    ! bonsai effective hit at lwatert*/
        float clwallsk;            /**** linfo=9: clwallsk*/
        float bswallsk;            /**** linfo=10: bswallsk*/
        float bsdir_lfdir2[3];     /**** linfo=11,12,13: bsdir_lfdir2*/
        float bsenergy_lfdir2;     /**** linfo=14: bsenergy_lfdir2*/
        float spa_random_pos[3];   /**** linfo=15,16,17: spa_random_pos*/
        float spa_random_loglike;  /**** linfo=18: spa_random_loglike*/
        float spa_random_resq;     /**** linfo=19: spa_random_resq*/
        float spa_random_dt;       /**** linfo=20: spa_random_dt*/
        float spa_random_dlt;      /**** linfo=21: spa_random_dlt*/
        float spa_random_dll;      /**** linfo=22: spa_random_dll*/
        float bsn20raw;            /**** linfo=23: bonsai N20raw*/
        float bsr02;               /**** linfo=24: bonsai R02*/
        float bsclik;              /**** linfo=25: bonsai clik*/
        float bsovaq;              /**** linfo=26: bonsai ovaq*/
        /**** linfo=27-33: Ariadne*/
        float adir[3];  // 27-29
        float amsg;     // 30
        float aratio;   // 31
        float anscat;   // 32
        float acosscat; // 33
        /**** linfo=34-46: solar analysis*/
        float poswal[3];       /* 34-36, poswal from effwal*/
        float idtgsk_wobadch;  /* 37, skdetsim output without badch mask*/
        float dir_slered[3];   /* 38-40, for reprocess/B8MC  direction*/
        float ene_slered;      /* 41,  for reprocess/B8MC  energy*/
        float effwal_slered;   /* 42,  for reprocess/B8MC  effective wall*/
        float ovaq_slered;     /* 43,  for reprocess/B8MC  ovaq*/
        int multi_spal;        /* 44,  for spallation by Scott Locke*/
        int wit_isgood;        /* 45,  for spallation by Scott Locke*/
        int cloud;             /* 46,  for spallation by Scott Locke*/
        /**** linfo=47-50: unused */
        int dummies1[4];
        /**** linfo=51-54: relic analysis*/
        float q50;    /* 51,  for relic analysis*/
        float angle;  /* 52,  for relic analysis*/
        float mult;   /* 53,  for relic analysis*/
        float pilike; /* 54,  for relic analysis*/
        /**** linfo 55-60: unused */
        int dummies2[6];
        /**** linfo=61-65: software trigger*/
        int swtrig;         /* 61, swtrig*/
        int swtrig_thr[4];  /* 62--65, threshold of trigid = 0--3 */
        /**** linfo=66--70: unused */
        int dummies3[5];
        /**** for T2K lowe (71-85)*/
        float t2kdt0;  // 71
        float t2kgpstvtx; // 72
        int t2kiseldt0; // 73
        int t2kn30max; // 74
        float t2kn30time; // 75
        /* reserved for the future */
        int t2kdummies[10]; // 76-85
        /**** save SK-I info in linfo(101-141)*/
        float kaienergy; // 101
        float kaivertex[3]; // 102-104
        float kaigoodness;  // 105
        float kaidir[3];    // 106-108
        float kaidirks;     // 109
        int kain30; // 110
        int kain50; // 111
        int kain200; // 112
        float kaipatlik; // 113
        int kaimaxnn; // 114
        int kainnall; // 115
        int kainume; // 116
        int kaindeno; // 117
        float clcl_vertex[4]; // 118-121
        float clcl_goodness; // 122
        float mrcvertex[3]; // 123-125
        float mrcgoodness; // 126
        /**** for muyn=0 spacut*/
        int spaevnum0; // 127
        float spaloglike0; // 128
        float sparesq0; // 129
        float spadt0; // 130
        float spadlt0; // 131
        int spamuyn0; // 132
        float spamugdn0; // 133
        /**** for B8MC*/
        int kailfflag; // 134
        int clcl_flag; // 135
        float kaieffhit[6]; // 136-141
        /**** linfo 142-256: unused */
        int dummies4[114];
    };
  };
} skroot_lowe_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct skroot_lowe_common *skroot_lowe;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct skroot_lowe_common *skroot_lowe = &skroot_lowe_;
#endif

extern struct skroot_mu_common {
  float  muentpoint[3];
  float  mudir[3];
  double mutimediff;
  float  mugoodness;
  float  muqismsk;
  int    muyn;
  int    mufast_flag;
  int    muboy_status;
  int    muboy_ntrack;
  float  muboy_entpos[10][4];
  float  muboy_dir[3];
  float  muboy_goodness;
  float  muboy_length;
  float  muboy_dedx[200];
  float  mubff_entpos[3];
  float  mubff_dir[3];
  float  mubff_goodness;
  int    muninfo;
  union {
      int muinfo[255];
      // muinfo=2: number of sub-trigger
      // muinfo=5: parent muon event number
      // muinfo=6: number of subtrigger in AFT
      float rmuinfo[255];
      // rmuinfo(1): original qismsk before fix_maxqisk()
      // rmuinfo(3): original qimxsk before fix_maxqisk()
      // and for some reason rmuinfo[4] is given it's own name
      struct {
          float muinfo_dummies_1[3];
          float muqimxsk;  // equivalent to muinfo[3], qimxsk after fix_maxqisk()
          int muinfo_dummies_2[251];
      };
  };
} skroot_mu_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct skroot_mu_common *skroot_mu;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct skroot_mu_common *skroot_mu = &skroot_mu_;
#endif

extern struct skroot_sle_common {
  float  itwallcut;
  int    itnsel;
  float  itbsvertex[4];  // bonsai fit
  float  itbsresult[6];
  float  itbsgood[3];
  int    itnbonsai;
  float  itcfgoodness;
  int    itnclusfit;
  /**** to keep SK-IV SLE reduction info.*/
  union {
      float itcfvertex[4];
      struct {
          float itbsdir[3];
          float itbsdirdummy;
      };
  };
  union {
      float itcfresult[4]; // clusfit
      struct {
          float itbsenergy;
          float itbsdirks;
          float itbswallsk;
          float itcfresultdummy;
      };
  };
} skroot_sle_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct skroot_sle_common *skroot_sle;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct skroot_sle_common *skroot_sle = &skroot_sle_;
#endif


/*------ data statements -----------------------------------------------*/


#ifndef NO_STATIC_DATA


#endif  /* #ifndef NO_STATIC_DATA */


/*------ end of fortran header -----------------------------------------*/


#ifdef __cplusplus
}
#endif


#endif  /* #ifndef FH2H_INC_SKROOT_LOWE_H */
