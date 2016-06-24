//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep  4 11:28:11 2015 by ROOT version 5.28/00c
// from TTree h1/DST
// found on file: jan14sk4_skdetsim13p90_neut532.reduc.081_fQv4r0.root
//////////////////////////////////////////////////////////

#ifndef fqEvent_h
#define fqEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

class fqEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nring;
   UInt_t          nrun;
   Int_t           nev;
   Int_t           nsub;
   UInt_t          cate;
   Float_t         potot;
   UInt_t          nhit;
   Float_t         pomax;
   Float_t         potota;
   UInt_t          nhita;
   UInt_t          nhitac;
   Float_t         pomaxa;
   Float_t         wall;
   Float_t         evis;
   Float_t         rtsum;
   Float_t         rtmax;
   Float_t         wlen;
   UInt_t          ip[10];   //[nring]
   Float_t         pos[3];
   Float_t         dir[10][3];   //[nring]
   Float_t         dirtot[3];
   Float_t         ang[10];   //[nring]
   Float_t         rtot[10];   //[nring]
   Float_t         amom[10];   //[nring]
   Float_t         rtote[10];   //[nring]
   Float_t         amome[10];   //[nring]
   Float_t         rtotm[10];   //[nring]
   Float_t         amomm[10];   //[nring]
   UInt_t          nsube;
   UInt_t          ndcy;
   UInt_t          ngate;
   UInt_t          nbye;
   Float_t         probms[10][6];   //[nring]
   Float_t         prmslg[10][6];   //[nring]
   Int_t           date[3];
   Int_t           time[4];
   Float_t         elpsday;
   Int_t           numpo[10];   //[nring]
   Float_t         apos[3];
   Float_t         adir[3];
   Float_t         aang;
   Float_t         agood;
   Float_t         wgain;
   Int_t           nbad;
   Int_t           nbada;
   Float_t         msdir[10][3][6];   //[nring]
   Float_t         amomp[10];
   Float_t         ange[10];
   Float_t         angm[10];
   Float_t         angp[10];
   Int_t           ntot[10];
   Float_t         probth[10][6];
   Float_t         probpt[10][6];
   Float_t         pi0like[2];
   Float_t         pi0_e[2][2];
   Float_t         pi0_dir[2][2][3];
   Float_t         pi0mass[2];
   Float_t         evisold;
   Float_t         evisoldxe;
   Float_t         evisnew;
   Int_t           nmue;
   UInt_t          etype[10];   //[nmue]
   Float_t         etime[10];   //[nmue]
   Float_t         epos[10][3];   //[nmue]
   Float_t         edir[10][3];   //[nmue]
   Float_t         egood[10];   //[nmue]
   Float_t         ehit[10];   //[nmue]
   Float_t         mueprob[2];
   Int_t           nmue_sel;
   Int_t           Rnring;
   Float_t         Rdir[30][3];   //[Rnring]
   Float_t         Rang[30];   //[Rnring]
   UInt_t          Riring;
   Float_t         Rtwout[30];   //[Rnring]
   Float_t         Rtwith[30];   //[Rnring]
   Float_t         Alwout;
   Float_t         Alwith;
   Float_t         Qsmi;
   Float_t         Qsmo;
   Float_t         Qexi;
   Float_t         Qexo;
   Float_t         Pe5d;
   Float_t         En5d;
   Float_t         Eh5d;
   Float_t         Pe5do;
   Float_t         En5do;
   Float_t         Eh5do;
   Float_t         Rtadd;
   Float_t         Pdgeta;
   Float_t         Pd5d;
   Float_t         Pdthre;
   Float_t         Pd5do;
   Float_t         Delpd;
   Float_t         Ropena[30];   //[Rnring]
   Int_t           Maxth;
   Float_t         Pkang;
   Float_t         Qrfct;
   Float_t         Pdfct;
   Float_t         Pkfct;
   Float_t         Agfct;
   Float_t         Dlfct;
   Int_t           Iflag;
   Float_t         Pmfct;
   Float_t         Imfct;
   Float_t         Rilike;
   Int_t           ri_ver;
   Float_t         ri_pid;
   Int_t           ri_nring;
   Float_t         ri_flag[10];   //[ri_nring]
   Float_t         ri_dlfct[10];   //[ri_nring]
   Float_t         ri_pdfct[10];   //[ri_nring]
   Float_t         ri_pkfct[10];   //[ri_nring]
   Float_t         ri_vafct[10];   //[ri_nring]
   Float_t         ri_total[10];   //[ri_nring]
   Float_t         ri_dir[10][3];   //[ri_nring]
   Float_t         ri_imfct[10];   //[ri_nring]
   Float_t         ri_pmfct[10];   //[ri_nring]
   Int_t           npar;
   Float_t         wallv;
   UInt_t          ipv[50];   //[npar]
   Float_t         posv[3];
   Float_t         dirv[50][3];   //[npar]
   Float_t         pmomv[50];   //[npar]
   Int_t           light_flag[50];   //[npar]
   Int_t           npar2;
   Float_t         wallv2[50];   //[npar2]
   UInt_t          ipv2[50];   //[npar2]
   UInt_t          iorg[50];   //[npar2]
   Float_t         posv2[50][3];   //[npar2]
   Float_t         dirv2[50][3];   //[npar2]
   Float_t         pmomv2[50];   //[npar2]
   Int_t           numnu;
   Int_t           mode;
   Int_t           ipnu[50];   //[numnu]
   Float_t         pnu[50];   //[numnu]
   Float_t         dirnu[50][3];   //[numnu]
   Float_t         flxg[3];
   Float_t         flxh01[3];
   Float_t         kflux[4];
   Float_t         bs71[3];
   Float_t         bs74[3];
   Float_t         flxf[3];
   Float_t         flxh1d[3];
   Float_t         flxb03[3];
   Float_t         flxf03[3];
   Float_t         flxh06[3];
   Float_t         flxh11[3];
   Int_t           scan[2];
   Int_t           nuceff_version;
   Int_t           charge_exchange;
   Int_t           absorbed;
   Int_t           multipi_gen;
   Int_t           scattering;
   Int_t           piless_dcy;
   Int_t           nscndprt;
   Int_t           itrkscnd[1000];   //[nscndprt]
   Int_t           istakscnd[1000];   //[nscndprt]
   Float_t         vtxscnd[1000][3];   //[nscndprt]
   Float_t         pscnd[1000][3];   //[nscndprt]
   Int_t           iprtscnd[1000];   //[nscndprt]
   Float_t         tscnd[1000];   //[nscndprt]
   Int_t           iprntprt[1000];   //[nscndprt]
   Int_t           lmecscnd[1000];   //[nscndprt]
   Int_t           iprnttrk[1000];   //[nscndprt]
   Int_t           iorgprt[1000];   //[nscndprt]
   Float_t         pprnt[1000][3];   //[nscndprt]
   Int_t           iflgscnd[1000];   //[nscndprt]
   Float_t         pprntinit[1000][3];   //[nscndprt]
   Float_t         vtxprnt[1000][3];   //[nscndprt]
   Int_t           iprntidx[1000];   //[nscndprt]
   Int_t           nchilds[1000];   //[nscndprt]
   Int_t           ichildidx[1000];   //[nscndprt]
   Int_t           Npvcscnd;
   Int_t           Nchildsvc[100];   //[Npvcscnd]
   Int_t           Ichildidxvc[100];   //[Npvcscnd]
   Int_t           Neutmode;
   Int_t           Npne;
   Int_t           Ipne[100];   //[Npne]
   Float_t         Abspne[100];   //[Npne]
   Float_t         Pvecne[100][3];   //[Npne]
   Int_t           Npvc;
   Int_t           Ipvc[100];   //[Npvc]
   Int_t           Ichvc[100];   //[Npvc]
   Int_t           Iorgvc[100];   //[Npvc]
   Int_t           Iflvc[100];   //[Npvc]
   Float_t         Abspvc[100];   //[Npvc]
   Float_t         Pvc[100][3];   //[Npvc]
   Float_t         Crsx;
   Float_t         Crsy;
   Float_t         Crsz;
   Float_t         Crsphi;
   Int_t           Numbndn;
   Int_t           Numbndp;
   Int_t           Numfrep;
   Int_t           Numatom;
   Int_t           Ibound;
   Int_t           Neutinpmt;
   Int_t           Nvert;
   Float_t         Posvert[300][3];   //[Nvert]
   Int_t           Iflgvert[300];   //[Nvert]
   Int_t           Nvcvert;
   Float_t         Dirvert[900][3];   //[Nvcvert]
   Float_t         Abspvert[900];   //[Nvcvert]
   Float_t         Abstpvert[900];   //[Nvcvert]
   Int_t           Ipvert[900];   //[Nvcvert]
   Int_t           Iverti[900];   //[Nvcvert]
   Int_t           Ivertf[900];   //[Nvcvert]
   Float_t         Fsiprob;
   Float_t         Fsivarwt[24];
   Int_t           Nvertsk;
   Float_t         Posvertsk[300][3];   //[Nvertsk]
   Int_t           Iflgvertsk[300];   //[Nvertsk]
   Int_t           Nvcvertsk;
   Float_t         Dirvertsk[900][3];   //[Nvcvertsk]
   Float_t         Abspvertsk[900];   //[Nvcvertsk]
   Float_t         Abstpvertsk[900];   //[Nvcvertsk]
   Int_t           Ipvertsk[900];   //[Nvcvertsk]
   Int_t           Ivertisk[900];   //[Nvcvertsk]
   Int_t           Ivertfsk[900];   //[Nvcvertsk]
   Float_t         Fsiprobsk;
   Float_t         trgoff;
   Int_t           fqntwnd;
   Int_t           fqtwnd_iclstr[10];   //[fqntwnd]
   Int_t           fqtwnd_npeak[10];   //[fqntwnd]
   Float_t         fqtwnd_prftt0[10];   //[fqntwnd]
   Float_t         fqtwnd_prftpos[10][3];   //[fqntwnd]
   Float_t         fqtwnd[10][2];   //[fqntwnd]
   Float_t         fqtwnd_peakt0[10][10];   //[fqntwnd]
   Float_t         fqtwnd_peakiness[10][10];   //[fqntwnd]
   Int_t           fqnse;
   Int_t           fqitwnd[10];   //[fqnse]
   Int_t           fqipeak[10];   //[fqnse]
   Int_t           fqnhitpmt[10];   //[fqnse]
   Float_t         fqtotq[10];   //[fqnse]
   Float_t         fq0rtotmu[10];   //[fqnse]
   Float_t         fq0rnll[10];   //[fqnse]
   Int_t           fqn50[10];   //[fqnse]
   Float_t         fqq50[10];   //[fqnse]
   Int_t           fq1rpcflg[10][7];   //[fqnse]
   Float_t         fq1rmom[10][7];   //[fqnse]
   Float_t         fq1rt0[10][7];   //[fqnse]
   Float_t         fq1rtotmu[10][7];   //[fqnse]
   Float_t         fq1rnll[10][7];   //[fqnse]
   Float_t         fq1rpos[10][7][3];   //[fqnse]
   Float_t         fq1rdir[10][7][3];   //[fqnse]
   Float_t         fq1rdconv[10][7];   //[fqnse]
   Float_t         fq1reloss[10][7];   //[fqnse]
   Int_t           fqpi0pcflg[2];
   Float_t         fqpi0mom1[2];
   Float_t         fqpi0mom2[2];
   Float_t         fqpi0momtot[2];
   Float_t         fqpi0dconv1[2];
   Float_t         fqpi0dconv2[2];
   Float_t         fqpi0t0[2];
   Float_t         fqpi0totmu[2];
   Float_t         fqpi0nll[2];
   Float_t         fqpi0mass[2];
   Float_t         fqpi0photangle[2];
   Float_t         fqpi0pos[2][3];
   Float_t         fqpi0dir1[2][3];
   Float_t         fqpi0dir2[2][3];
   Float_t         fqpi0dirtot[2][3];
   Int_t           fqnmrfit;
   Int_t           fqmrifit[32];   //[fqnmrfit]
   Int_t           fqmrnring[32];   //[fqnmrfit]
   Int_t           fqmrpcflg[32];   //[fqnmrfit]
   Float_t         fqmrnll[32];   //[fqnmrfit]
   Float_t         fqmrtotmu[32];   //[fqnmrfit]
   Int_t           fqmrpid[32][6];   //[fqnmrfit]
   Float_t         fqmrmom[32][6];   //[fqnmrfit]
   Float_t         fqmrdconv[32][6];   //[fqnmrfit]
   Float_t         fqmreloss[32][6];   //[fqnmrfit]
   Float_t         fqmrt0[32][6];   //[fqnmrfit]
   Float_t         fqmrpos[32][6][3];   //[fqnmrfit]
   Float_t         fqmrdir[32][6][3];   //[fqnmrfit]
   Int_t           fqmsnfit;
   Int_t           fqmspcflg[5];   //[fqmsnfit]
   Int_t           fqmsnseg[5];   //[fqmsnfit]
   Int_t           fqmspid[5];   //[fqmsnfit]
   Int_t           fqmsifit[5];   //[fqmsnfit]
   Int_t           fqmsimer[5];   //[fqmsnfit]
   Float_t         fqmstotmu[5];   //[fqmsnfit]
   Float_t         fqmsnll[5];   //[fqmsnfit]
   Float_t         fqmsmom[5][20];   //[fqmsnfit]
   Float_t         fqmseloss[5][20];   //[fqmsnfit]
   Float_t         fqmst0[5][20];   //[fqmsnfit]
   Float_t         fqmspos[5][20][3];   //[fqmsnfit]
   Float_t         fqmsdir[5][20][3];   //[fqmsnfit]
   Float_t         prid_ang;
   Float_t         prid_len;
   Float_t         prid_pmom;
   Float_t         prid_mang;
   Float_t         prid_plike;
   Float_t         prid_mlike;
   Int_t           prid_nused;
   Int_t           nn;
   Int_t           trgtype;
   Int_t           n200m;
   Float_t         t200m;
   Float_t         lasthit;
   Float_t         ntime[1000];   //[nn]
   Float_t         ngood[1000];   //[nn]
   Float_t         nvx[1000][3];   //[nn]
   Float_t         bvx[1000][3];   //[nn]
   Float_t         dirtotepi[3];
   Float_t         dirtotenpi[3];
   Float_t         dirtotmue[3];
   Float_t         dirsum[3];
   Float_t         etot;
   Float_t         etotepi;
   Float_t         etotenpi;
   Float_t         etotmue;
   Float_t         oscweight[2][4];
   Float_t         oscwgt;
   Float_t         live;
   Float_t         sacth;
   Float_t         sactg;
   Float_t         sacth1d;
   Float_t         ccqelike;
   Float_t         mulike;
   Float_t         mproton;
   Float_t         outdir[3];
   Int_t           nused2;
   Float_t         vsum;
   Float_t         qsep;
   Float_t         qpeak;
   Float_t         qhigh;
   Float_t         qlow;
   Int_t           nzeros2;
   Float_t         moyenne2;
   Float_t         lproton;
   Float_t         qwin;
   Float_t         ent_pos[3];
   Float_t         ent_dir[3];
   Float_t         length;
   Float_t         tr_mom1;
   Float_t         A_ent_mom;
   Float_t         A_ent_pos[3];
   Float_t         A_ent_dir[3];
   Float_t         A_ext_mom;
   Float_t         A_ext_pos[3];
   Float_t         A_ext_dir[3];
   Float_t         Fit_pos[3];
   Float_t         Fit_dir[3];
   Float_t         Fit_len;
   Float_t         Fit_mom;
   Int_t           Fit_pid;
   Int_t           Um_ehit8m;
   Int_t           Um_ohit8m;
   Float_t         Um_qent;
   Float_t         Sh_chi1p;
   Float_t         Sh_delta;
   Float_t         Sh_mean;
   Float_t         Sh_meanq;
   Int_t           Sh_id;
   Float_t         Tr_stop[3];
   Float_t         Tr_mom;
   Float_t         Tr_len;
   Float_t         Tr_len1;
   Int_t           Pid_flg;
   Float_t         Crs1;
   Float_t         Crs2;
   Int_t           iclass;
   Int_t           mu_class;
   Int_t           mu_dec;
   Float_t         mu_dir[3];
   Float_t         mu_pos[3];
   Float_t         mu_good;
   Int_t           history;
   Int_t           Pdst;
   Int_t           idoff;
   Float_t         anthit;
   Int_t           idseq;
   Float_t         tstfrac;
   Int_t           judge;
   Float_t         Upcrs1;
   Float_t         Upcrs2;
   Float_t         lst;
   Int_t           jd;
   Float_t         fjd;
   Float_t         alt;
   Float_t         azi;
   Float_t         ra;
   Float_t         dec;
   Float_t         glat;
   Float_t         glong;
   Int_t           t2knuefs;
   Int_t           t2knuenpi0;
   Int_t           t2knueflag;
   Double_t        totwgt;



   // List of branches
   TBranch        *b_nring;   //!
   TBranch        *b_nrun;   //!
   TBranch        *b_nev;   //!
   TBranch        *b_nsub;   //!
   TBranch        *b_cate;   //!
   TBranch        *b_potot;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_pomax;   //!
   TBranch        *b_potota;   //!
   TBranch        *b_nhita;   //!
   TBranch        *b_nhitac;   //!
   TBranch        *b_pomaxa;   //!
   TBranch        *b_wall;   //!
   TBranch        *b_evis;   //!
   TBranch        *b_rtsum;   //!
   TBranch        *b_rtmax;   //!
   TBranch        *b_wlen;   //!
   TBranch        *b_ip;   //!
   TBranch        *b_pos;   //!
   TBranch        *b_dir;   //!
   TBranch        *b_dirtot;   //!
   TBranch        *b_ang;   //!
   TBranch        *b_rtot;   //!
   TBranch        *b_amom;   //!
   TBranch        *b_rtote;   //!
   TBranch        *b_amome;   //!
   TBranch        *b_rtotm;   //!
   TBranch        *b_amomm;   //!
   TBranch        *b_nsube;   //!
   TBranch        *b_ndcy;   //!
   TBranch        *b_ngate;   //!
   TBranch        *b_nbye;   //!
   TBranch        *b_probms;   //!
   TBranch        *b_prmslg;   //!
   TBranch        *b_date;   //!
   TBranch        *b_time;   //!
   TBranch        *b_elpsday;   //!
   TBranch        *b_numpo;   //!
   TBranch        *b_apos;   //!
   TBranch        *b_adir;   //!
   TBranch        *b_aang;   //!
   TBranch        *b_agood;   //!
   TBranch        *b_wgain;   //!
   TBranch        *b_nbad;   //!
   TBranch        *b_nbada;   //!
   TBranch        *b_msdir;   //!
   TBranch        *b_amomp;   //!
   TBranch        *b_ange;   //!
   TBranch        *b_angm;   //!
   TBranch        *b_angp;   //!
   TBranch        *b_ntot;   //!
   TBranch        *b_probth;   //!
   TBranch        *b_probpt;   //!
   TBranch        *b_pi0like;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_dir;   //!
   TBranch        *b_pi0mass;   //!
   TBranch        *b_evisold;   //!
   TBranch        *b_evisoldxe;   //!
   TBranch        *b_evisnew;   //!
   TBranch        *b_nmue;   //!
   TBranch        *b_etype;   //!
   TBranch        *b_etime;   //!
   TBranch        *b_epos;   //!
   TBranch        *b_edir;   //!
   TBranch        *b_egood;   //!
   TBranch        *b_ehit;   //!
   TBranch        *b_mueprob;   //!
   TBranch        *b_nmue_sel;   //!
   TBranch        *b_Rnring;   //!
   TBranch        *b_Rdir;   //!
   TBranch        *b_Rang;   //!
   TBranch        *b_Riring;   //!
   TBranch        *b_Rtwout;   //!
   TBranch        *b_Rtwith;   //!
   TBranch        *b_Alwout;   //!
   TBranch        *b_Alwith;   //!
   TBranch        *b_Qsmi;   //!
   TBranch        *b_Qsmo;   //!
   TBranch        *b_Qexi;   //!
   TBranch        *b_Qexo;   //!
   TBranch        *b_Pe5d;   //!
   TBranch        *b_En5d;   //!
   TBranch        *b_Eh5d;   //!
   TBranch        *b_Pe5do;   //!
   TBranch        *b_En5do;   //!
   TBranch        *b_Eh5do;   //!
   TBranch        *b_Rtadd;   //!
   TBranch        *b_Pdgeta;   //!
   TBranch        *b_Pd5d;   //!
   TBranch        *b_Pdthre;   //!
   TBranch        *b_Pd5do;   //!
   TBranch        *b_Delpd;   //!
   TBranch        *b_Ropena;   //!
   TBranch        *b_Maxth;   //!
   TBranch        *b_Pkang;   //!
   TBranch        *b_Qrfct;   //!
   TBranch        *b_Pdfct;   //!
   TBranch        *b_Pkfct;   //!
   TBranch        *b_Agfct;   //!
   TBranch        *b_Dlfct;   //!
   TBranch        *b_Iflag;   //!
   TBranch        *b_Pmfct;   //!
   TBranch        *b_Imfct;   //!
   TBranch        *b_Rilike;   //!
   TBranch        *b_ri_ver;   //!
   TBranch        *b_ri_pid;   //!
   TBranch        *b_ri_nring;   //!
   TBranch        *b_ri_flag;   //!
   TBranch        *b_ri_dlfct;   //!
   TBranch        *b_ri_pdfct;   //!
   TBranch        *b_ri_pkfct;   //!
   TBranch        *b_ri_vafct;   //!
   TBranch        *b_ri_total;   //!
   TBranch        *b_ri_dir;   //!
   TBranch        *b_ri_imfct;   //!
   TBranch        *b_ri_pmfct;   //!
   TBranch        *b_npar;   //!
   TBranch        *b_wallv;   //!
   TBranch        *b_ipv;   //!
   TBranch        *b_posv;   //!
   TBranch        *b_dirv;   //!
   TBranch        *b_pmomv;   //!
   TBranch        *b_light_flag;   //!
   TBranch        *b_npar2;   //!
   TBranch        *b_wallv2;   //!
   TBranch        *b_ipv2;   //!
   TBranch        *b_iorg;   //!
   TBranch        *b_posv2;   //!
   TBranch        *b_dirv2;   //!
   TBranch        *b_pmomv2;   //!
   TBranch        *b_numnu;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_ipnu;   //!
   TBranch        *b_pnu;   //!
   TBranch        *b_dirnu;   //!
   TBranch        *b_flxg;   //!
   TBranch        *b_flxh01;   //!
   TBranch        *b_kflux;   //!
   TBranch        *b_bs71;   //!
   TBranch        *b_bs74;   //!
   TBranch        *b_flxf;   //!
   TBranch        *b_flxh1d;   //!
   TBranch        *b_flxb03;   //!
   TBranch        *b_flxf03;   //!
   TBranch        *b_flxh06;   //!
   TBranch        *b_flxh11;   //!
   TBranch        *b_scan;   //!
   TBranch        *b_nuceff_version;   //!
   TBranch        *b_charge_exchange;   //!
   TBranch        *b_absorbed;   //!
   TBranch        *b_multipi_gen;   //!
   TBranch        *b_scattering;   //!
   TBranch        *b_piless_dcy;   //!
   TBranch        *b_nscndprt;   //!
   TBranch        *b_itrkscnd;   //!
   TBranch        *b_istakscnd;   //!
   TBranch        *b_vtxscnd;   //!
   TBranch        *b_pscnd;   //!
   TBranch        *b_iprtscnd;   //!
   TBranch        *b_tscnd;   //!
   TBranch        *b_iprntprt;   //!
   TBranch        *b_lmecscnd;   //!
   TBranch        *b_iprnttrk;   //!
   TBranch        *b_iorgprt;   //!
   TBranch        *b_pprnt;   //!
   TBranch        *b_iflgscnd;   //!
   TBranch        *b_pprntinit;   //!
   TBranch        *b_vtxprnt;   //!
   TBranch        *b_iprntidx;   //!
   TBranch        *b_nchilds;   //!
   TBranch        *b_ichildidx;   //!
   TBranch        *b_Npvcscnd;   //!
   TBranch        *b_Nchildsvc;   //!
   TBranch        *b_Ichildidxvc;   //!
   TBranch        *b_Neutmode;   //!
   TBranch        *b_Npne;   //!
   TBranch        *b_Ipne;   //!
   TBranch        *b_Abspne;   //!
   TBranch        *b_Pvecne;   //!
   TBranch        *b_Npvc;   //!
   TBranch        *b_Ipvc;   //!
   TBranch        *b_Ichvc;   //!
   TBranch        *b_Iorgvc;   //!
   TBranch        *b_Iflvc;   //!
   TBranch        *b_Abspvc;   //!
   TBranch        *b_Pvc;   //!
   TBranch        *b_Crsx;   //!
   TBranch        *b_Crsy;   //!
   TBranch        *b_Crsz;   //!
   TBranch        *b_Crsphi;   //!
   TBranch        *b_Numbndn;   //!
   TBranch        *b_Numbndp;   //!
   TBranch        *b_Numfrep;   //!
   TBranch        *b_Numatom;   //!
   TBranch        *b_Ibound;   //!
   TBranch        *b_Neutinpmt;   //!
   TBranch        *b_Nvert;   //!
   TBranch        *b_Posvert;   //!
   TBranch        *b_Iflgvert;   //!
   TBranch        *b_Nvcvert;   //!
   TBranch        *b_Dirvert;   //!
   TBranch        *b_Abspvert;   //!
   TBranch        *b_Abstpvert;   //!
   TBranch        *b_Ipvert;   //!
   TBranch        *b_Iverti;   //!
   TBranch        *b_Ivertf;   //!
   TBranch        *b_Fsiprob;   //!
   TBranch        *b_Fsivarwt;   //!
   TBranch        *b_Nvertsk;   //!
   TBranch        *b_Posvertsk;   //!
   TBranch        *b_Iflgvertsk;   //!
   TBranch        *b_Nvcvertsk;   //!
   TBranch        *b_Dirvertsk;   //!
   TBranch        *b_Abspvertsk;   //!
   TBranch        *b_Abstpvertsk;   //!
   TBranch        *b_Ipvertsk;   //!
   TBranch        *b_Ivertisk;   //!
   TBranch        *b_Ivertfsk;   //!
   TBranch        *b_Fsiprobsk;   //!
   TBranch        *b_trgoff;   //!
   TBranch        *b_fqntwnd;   //!
   TBranch        *b_fqtwnd_iclstr;   //!
   TBranch        *b_fqtwnd_npeak;   //!
   TBranch        *b_fqtwnd_prftt0;   //!
   TBranch        *b_fqtwnd_prftpos;   //!
   TBranch        *b_fqtwnd;   //!
   TBranch        *b_fqtwnd_peakt0;   //!
   TBranch        *b_fqtwnd_peakiness;   //!
   TBranch        *b_fqnse;   //!
   TBranch        *b_fqitwnd;   //!
   TBranch        *b_fqipeak;   //!
   TBranch        *b_fqnhitpmt;   //!
   TBranch        *b_fqtotq;   //!
   TBranch        *b_fq0rtotmu;   //!
   TBranch        *b_fq0rnll;   //!
   TBranch        *b_fqn50;   //!
   TBranch        *b_fqq50;   //!
   TBranch        *b_fq1rpcflg;   //!
   TBranch        *b_fq1rmom;   //!
   TBranch        *b_fq1rt0;   //!
   TBranch        *b_fq1rtotmu;   //!
   TBranch        *b_fq1rnll;   //!
   TBranch        *b_fq1rpos;   //!
   TBranch        *b_fq1rdir;   //!
   TBranch        *b_fq1rdconv;   //!
   TBranch        *b_fq1reloss;   //!
   TBranch        *b_fqpi0pcflg;   //!
   TBranch        *b_fqpi0mom1;   //!
   TBranch        *b_fqpi0mom2;   //!
   TBranch        *b_fqpi0momtot;   //!
   TBranch        *b_fqpi0dconv1;   //!
   TBranch        *b_fqpi0dconv2;   //!
   TBranch        *b_fqpi0t0;   //!
   TBranch        *b_fqpi0totmu;   //!
   TBranch        *b_fqpi0nll;   //!
   TBranch        *b_fqpi0mass;   //!
   TBranch        *b_fqpi0photangle;   //!
   TBranch        *b_fqpi0pos;   //!
   TBranch        *b_fqpi0dir1;   //!
   TBranch        *b_fqpi0dir2;   //!
   TBranch        *b_fqpi0dirtot;   //!
   TBranch        *b_fqnmrfit;   //!
   TBranch        *b_fqmrifit;   //!
   TBranch        *b_fqmrnring;   //!
   TBranch        *b_fqmrpcflg;   //!
   TBranch        *b_fqmrnll;   //!
   TBranch        *b_fqmrtotmu;   //!
   TBranch        *b_fqmrpid;   //!
   TBranch        *b_fqmrmom;   //!
   TBranch        *b_fqmrdconv;   //!
   TBranch        *b_fqmreloss;   //!
   TBranch        *b_fqmrt0;   //!
   TBranch        *b_fqmrpos;   //!
   TBranch        *b_fqmrdir;   //!
   TBranch        *b_fqmsnfit;   //!
   TBranch        *b_fqmspcflg;   //!
   TBranch        *b_fqmsnseg;   //!
   TBranch        *b_fqmspid;   //!
   TBranch        *b_fqmsifit;   //!
   TBranch        *b_fqmsimer;   //!
   TBranch        *b_fqmstotmu;   //!
   TBranch        *b_fqmsnll;   //!
   TBranch        *b_fqmsmom;   //!
   TBranch        *b_fqmseloss;   //!
   TBranch        *b_fqmst0;   //!
   TBranch        *b_fqmspos;   //!
   TBranch        *b_fqmsdir;   //!
   TBranch        *b_prid_ang;   //!
   TBranch        *b_prid_len;   //!
   TBranch        *b_prid_pmom;   //!
   TBranch        *b_prid_mang;   //!
   TBranch        *b_prid_plike;   //!
   TBranch        *b_prid_mlike;   //!
   TBranch        *b_prid_nused;   //!
   TBranch        *b_nn;   //!
   TBranch        *b_trgtype;   //!
   TBranch        *b_n200m;   //!
   TBranch        *b_t200m;   //!
   TBranch        *b_lasthit;   //!
   TBranch        *b_ntime;   //!
   TBranch        *b_ngood;   //!
   TBranch        *b_nvx;   //!
   TBranch        *b_bvx;   //!
   TBranch        *b_dirtotepi;   //!
   TBranch        *b_dirtotenpi;   //!
   TBranch        *b_dirtotmue;   //!
   TBranch        *b_dirsum;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_etotepi;   //!
   TBranch        *b_etotenpi;   //!
   TBranch        *b_etotmue;   //!
   TBranch        *b_oscweight;   //!
   TBranch        *b_oscwgt;   //!
   TBranch        *b_live;   //!
   TBranch        *b_sacth;   //!
   TBranch        *b_sactg;   //!
   TBranch        *b_sacth1d;   //!
   TBranch        *b_ccqelike;   //!
   TBranch        *b_mulike;   //!
   TBranch        *b_mproton;   //!
   TBranch        *b_outdir;   //!
   TBranch        *b_nused2;   //!
   TBranch        *b_vsum;   //!
   TBranch        *b_qsep;   //!
   TBranch        *b_qpeak;   //!
   TBranch        *b_qhigh;   //!
   TBranch        *b_qlow;   //!
   TBranch        *b_nzeros2;   //!
   TBranch        *b_moyenne2;   //!
   TBranch        *b_lproton;   //!
   TBranch        *b_qwin;   //!
   TBranch        *b_ent_pos;   //!
   TBranch        *b_ent_dir;   //!
   TBranch        *b_length;   //!
   TBranch        *b_tr_mom1;   //!
   TBranch        *b_A_ent_mom;   //!
   TBranch        *b_A_ent_pos;   //!
   TBranch        *b_A_ent_dir;   //!
   TBranch        *b_A_ext_mom;   //!
   TBranch        *b_A_ext_pos;   //!
   TBranch        *b_A_ext_dir;   //!
   TBranch        *b_Fit_pos;   //!
   TBranch        *b_Fit_dir;   //!
   TBranch        *b_Fit_len;   //!
   TBranch        *b_Fit_mom;   //!
   TBranch        *b_Fit_pid;   //!
   TBranch        *b_Um_ehit8m;   //!
   TBranch        *b_Um_ohit8m;   //!
   TBranch        *b_Um_qent;   //!
   TBranch        *b_Sh_chi1p;   //!
   TBranch        *b_Sh_delta;   //!
   TBranch        *b_Sh_mean;   //!
   TBranch        *b_Sh_meanq;   //!
   TBranch        *b_Sh_id;   //!
   TBranch        *b_Tr_stop;   //!
   TBranch        *b_Tr_mom;   //!
   TBranch        *b_Tr_len;   //!
   TBranch        *b_Tr_len1;   //!
   TBranch        *b_Pid_flg;   //!
   TBranch        *b_Crs1;   //!
   TBranch        *b_Crs2;   //!
   TBranch        *b_iclass;   //!
   TBranch        *b_mu_class;   //!
   TBranch        *b_mu_dec;   //!
   TBranch        *b_mu_dir;   //!
   TBranch        *b_mu_pos;   //!
   TBranch        *b_mu_good;   //!
   TBranch        *b_history;   //!
   TBranch        *b_Pdst;   //!
   TBranch        *b_idoff;   //!
   TBranch        *b_anthit;   //!
   TBranch        *b_idseq;   //!
   TBranch        *b_tstfrac;   //!
   TBranch        *b_judge;   //!
   TBranch        *b_Upcrs1;   //!
   TBranch        *b_Upcrs2;   //!
   TBranch        *b_lst;   //!
   TBranch        *b_jd;   //!
   TBranch        *b_fjd;   //!
   TBranch        *b_alt;   //!
   TBranch        *b_azi;   //!
   TBranch        *b_ra;   //!
   TBranch        *b_dec;   //!
   TBranch        *b_glat;   //!
   TBranch        *b_glong;   //!
   TBranch        *b_t2knuefs;   //!
   TBranch        *b_t2knuenpi0;   //!
   TBranch        *b_t2knueflag;   //!
   TBranch        *b_totwgt;   //!

   fqEvent(TTree *tree=0,const char* ntuple_type = "");
   virtual ~fqEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree,const char* ntuple_type = "");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};


fqEvent::fqEvent(TTree *tree, const char* ntupletype)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("jan14sk4_skdetsim13p90_neut532.reduc.081_fQv4r0.root");
      if (!f) {
         f = new TFile("jan14sk4_skdetsim13p90_neut532.reduc.081_fQv4r0.root");
      }
      tree = (TTree*)gDirectory->Get("h1");

   }
   Init(tree);
}

fqEvent::~fqEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fqEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fqEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fqEvent::Init(TTree *tree, const char* ntuple_type)
{

   // specify the type of the ntuple to set addresses for
   TString ntype = ntuple_type; 
   if (!ntype.CompareTo("")){ ntype = "Atmospheric"; }

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.;
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


   // for Atmoserphic data and MC
   if (!ntype.CompareTo("Atmospheric")){
   fChain->SetBranchAddress("nring", &nring, &b_nring);
   fChain->SetBranchAddress("nrun", &nrun, &b_nrun);
   fChain->SetBranchAddress("nev", &nev, &b_nev);
   fChain->SetBranchAddress("nsub", &nsub, &b_nsub);
   fChain->SetBranchAddress("potot", &potot, &b_potot);
   fChain->SetBranchAddress("nhit", &nhit, &b_nhit);
   fChain->SetBranchAddress("pomax", &pomax, &b_pomax);
   fChain->SetBranchAddress("nhita", &nhita, &b_nhita);
   fChain->SetBranchAddress("nhitac", &nhitac, &b_nhitac);
   fChain->SetBranchAddress("wall", &wall, &b_wall);
   fChain->SetBranchAddress("evis", &evis, &b_evis);
   fChain->SetBranchAddress("wlen", &wlen, &b_wlen);
   fChain->SetBranchAddress("ip", ip, &b_ip);
   fChain->SetBranchAddress("pos", pos, &b_pos);
   fChain->SetBranchAddress("dir", dir, &b_dir);
   fChain->SetBranchAddress("dirtot", dirtot, &b_dirtot);
   fChain->SetBranchAddress("wgain", &wgain, &b_wgain);
   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
   fChain->SetBranchAddress("ipv", ipv, &b_ipv);
   fChain->SetBranchAddress("posv", posv, &b_posv);
   fChain->SetBranchAddress("dirv", dirv, &b_dirv);
   fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
   fChain->SetBranchAddress("light_flag", light_flag, &b_light_flag);
   fChain->SetBranchAddress("npar2", &npar2, &b_npar2);
   fChain->SetBranchAddress("wallv2", wallv2, &b_wallv2);
   fChain->SetBranchAddress("ipv2", ipv2, &b_ipv2);
   fChain->SetBranchAddress("iorg", iorg, &b_iorg);
   fChain->SetBranchAddress("posv2", posv2, &b_posv2);
   fChain->SetBranchAddress("dirv2", dirv2, &b_dirv2);
   fChain->SetBranchAddress("pmomv2", pmomv2, &b_pmomv2);
   fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("ipnu", ipnu, &b_ipnu);
   fChain->SetBranchAddress("pnu", pnu, &b_pnu);
   fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
   fChain->SetBranchAddress("fqntwnd", &fqntwnd, &b_fqntwnd);
   fChain->SetBranchAddress("fqtwnd_iclstr", fqtwnd_iclstr, &b_fqtwnd_iclstr);
   fChain->SetBranchAddress("fqtwnd_npeak", fqtwnd_npeak, &b_fqtwnd_npeak);
   fChain->SetBranchAddress("fqtwnd_prftt0", fqtwnd_prftt0, &b_fqtwnd_prftt0);
   fChain->SetBranchAddress("fqtwnd_prftpos", fqtwnd_prftpos, &b_fqtwnd_prftpos);
   fChain->SetBranchAddress("fqtwnd", fqtwnd, &b_fqtwnd);
   fChain->SetBranchAddress("fqtwnd_peakt0", fqtwnd_peakt0, &b_fqtwnd_peakt0);
   fChain->SetBranchAddress("fqtwnd_peakiness", fqtwnd_peakiness, &b_fqtwnd_peakiness);
   fChain->SetBranchAddress("fqnse", &fqnse, &b_fqnse);
   fChain->SetBranchAddress("fqitwnd", fqitwnd, &b_fqitwnd);
   fChain->SetBranchAddress("fqipeak", fqipeak, &b_fqipeak);
   fChain->SetBranchAddress("fqnhitpmt", fqnhitpmt, &b_fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", fqtotq, &b_fqtotq);
   fChain->SetBranchAddress("fq0rtotmu", fq0rtotmu, &b_fq0rtotmu);
   fChain->SetBranchAddress("fq0rnll", fq0rnll, &b_fq0rnll);
   fChain->SetBranchAddress("fqn50", fqn50, &b_fqn50);
   fChain->SetBranchAddress("fqq50", fqq50, &b_fqq50);
   fChain->SetBranchAddress("fq1rpcflg", fq1rpcflg, &b_fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", fq1rmom, &b_fq1rmom);
   fChain->SetBranchAddress("fq1rt0", fq1rt0, &b_fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", fq1rtotmu, &b_fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", fq1rnll, &b_fq1rnll);
   fChain->SetBranchAddress("fq1rpos", fq1rpos, &b_fq1rpos);
   fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
   fChain->SetBranchAddress("fq1rdconv", fq1rdconv, &b_fq1rdconv);
   fChain->SetBranchAddress("fq1reloss", fq1reloss, &b_fq1reloss);
   fChain->SetBranchAddress("fqpi0pcflg", fqpi0pcflg, &b_fqpi0pcflg);
   fChain->SetBranchAddress("fqpi0mom1", fqpi0mom1, &b_fqpi0mom1);
   fChain->SetBranchAddress("fqpi0mom2", fqpi0mom2, &b_fqpi0mom2);
   fChain->SetBranchAddress("fqpi0momtot", fqpi0momtot, &b_fqpi0momtot);
   fChain->SetBranchAddress("fqpi0dconv1", fqpi0dconv1, &b_fqpi0dconv1);
   fChain->SetBranchAddress("fqpi0dconv2", fqpi0dconv2, &b_fqpi0dconv2);
   fChain->SetBranchAddress("fqpi0t0", fqpi0t0, &b_fqpi0t0);
   fChain->SetBranchAddress("fqpi0totmu", fqpi0totmu, &b_fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", fqpi0nll, &b_fqpi0nll);
   fChain->SetBranchAddress("fqpi0mass", fqpi0mass, &b_fqpi0mass);
   fChain->SetBranchAddress("fqpi0photangle", fqpi0photangle, &b_fqpi0photangle);
   fChain->SetBranchAddress("fqpi0pos", fqpi0pos, &b_fqpi0pos);
   fChain->SetBranchAddress("fqpi0dir1", fqpi0dir1, &b_fqpi0dir1);
   fChain->SetBranchAddress("fqpi0dir2", fqpi0dir2, &b_fqpi0dir2);
   fChain->SetBranchAddress("fqpi0dirtot", fqpi0dirtot, &b_fqpi0dirtot);
   fChain->SetBranchAddress("fqnmrfit", &fqnmrfit, &b_fqnmrfit);
   fChain->SetBranchAddress("fqmrifit", fqmrifit, &b_fqmrifit);
   fChain->SetBranchAddress("fqmrnring", fqmrnring, &b_fqmrnring);
   fChain->SetBranchAddress("fqmrpcflg", fqmrpcflg, &b_fqmrpcflg);
   fChain->SetBranchAddress("fqmrnll", fqmrnll, &b_fqmrnll);
   fChain->SetBranchAddress("fqmrtotmu", fqmrtotmu, &b_fqmrtotmu);
   fChain->SetBranchAddress("fqmrpid", fqmrpid, &b_fqmrpid);
   fChain->SetBranchAddress("fqmrmom", fqmrmom, &b_fqmrmom);
   fChain->SetBranchAddress("fqmrdconv", fqmrdconv, &b_fqmrdconv);
   fChain->SetBranchAddress("fqmreloss", fqmreloss, &b_fqmreloss);
   fChain->SetBranchAddress("fqmrt0", fqmrt0, &b_fqmrt0);
   fChain->SetBranchAddress("fqmrpos", fqmrpos, &b_fqmrpos);
   fChain->SetBranchAddress("fqmrdir", fqmrdir, &b_fqmrdir);
   fChain->SetBranchAddress("fqmsnfit", &fqmsnfit, &b_fqmsnfit);
   fChain->SetBranchAddress("fqmspcflg", fqmspcflg, &b_fqmspcflg);
   fChain->SetBranchAddress("fqmsnseg", fqmsnseg, &b_fqmsnseg);
   fChain->SetBranchAddress("fqmspid", fqmspid, &b_fqmspid);
   fChain->SetBranchAddress("fqmsifit", fqmsifit, &b_fqmsifit);
   fChain->SetBranchAddress("fqmsimer", fqmsimer, &b_fqmsimer);
   fChain->SetBranchAddress("fqmstotmu", fqmstotmu, &b_fqmstotmu);
   fChain->SetBranchAddress("fqmsnll", fqmsnll, &b_fqmsnll);
   fChain->SetBranchAddress("fqmsmom", fqmsmom, &b_fqmsmom);
   fChain->SetBranchAddress("fqmseloss", fqmseloss, &b_fqmseloss);
   fChain->SetBranchAddress("fqmst0", fqmst0, &b_fqmst0);
   fChain->SetBranchAddress("fqmspos", fqmspos, &b_fqmspos);
   fChain->SetBranchAddress("fqmsdir", fqmsdir, &b_fqmsdir);
   }
   if (!ntype.CompareTo("T2KMCReduced")){
   fChain->SetBranchAddress("nring", &nring, &b_nring);
   fChain->SetBranchAddress("nev", &nev, &b_nev);
   fChain->SetBranchAddress("nhitac", &nhitac, &b_nhitac);
   fChain->SetBranchAddress("numpo", numpo, &b_numpo);
   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
   fChain->SetBranchAddress("ipv", ipv, &b_ipv);
   fChain->SetBranchAddress("posv", posv, &b_posv);
   fChain->SetBranchAddress("dirv", dirv, &b_dirv);
   fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
   fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("ipnu", ipnu, &b_ipnu);
   fChain->SetBranchAddress("pnu", pnu, &b_pnu);
   fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
   fChain->SetBranchAddress("nuceff_version", &nuceff_version, &b_nuceff_version);
   fChain->SetBranchAddress("fqntwnd", &fqntwnd, &b_fqntwnd);
   fChain->SetBranchAddress("fqtwnd_iclstr", fqtwnd_iclstr, &b_fqtwnd_iclstr);
   fChain->SetBranchAddress("fqtwnd_npeak", fqtwnd_npeak, &b_fqtwnd_npeak);
   fChain->SetBranchAddress("fqtwnd_prftt0", fqtwnd_prftt0, &b_fqtwnd_prftt0);
   fChain->SetBranchAddress("fqtwnd_prftpos", fqtwnd_prftpos, &b_fqtwnd_prftpos);
   fChain->SetBranchAddress("fqtwnd", fqtwnd, &b_fqtwnd);
   fChain->SetBranchAddress("fqtwnd_peakt0", fqtwnd_peakt0, &b_fqtwnd_peakt0);
   fChain->SetBranchAddress("fqtwnd_peakiness", fqtwnd_peakiness, &b_fqtwnd_peakiness);
   fChain->SetBranchAddress("fqnse", &fqnse, &b_fqnse);
   fChain->SetBranchAddress("fqitwnd", fqitwnd, &b_fqitwnd);
   fChain->SetBranchAddress("fqipeak", fqipeak, &b_fqipeak);
   fChain->SetBranchAddress("fqnhitpmt", fqnhitpmt, &b_fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", fqtotq, &b_fqtotq);
   fChain->SetBranchAddress("fq0rtotmu", fq0rtotmu, &b_fq0rtotmu);
   fChain->SetBranchAddress("fq0rnll", fq0rnll, &b_fq0rnll);
   fChain->SetBranchAddress("fqn50", fqn50, &b_fqn50);
   fChain->SetBranchAddress("fqq50", fqq50, &b_fqq50);
   fChain->SetBranchAddress("fq1rpcflg", fq1rpcflg, &b_fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", fq1rmom, &b_fq1rmom);
   fChain->SetBranchAddress("fq1rt0", fq1rt0, &b_fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", fq1rtotmu, &b_fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", fq1rnll, &b_fq1rnll);
   fChain->SetBranchAddress("fq1rpos", fq1rpos, &b_fq1rpos);
   fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
   fChain->SetBranchAddress("fq1rdconv", fq1rdconv, &b_fq1rdconv);
   fChain->SetBranchAddress("fq1reloss", fq1reloss, &b_fq1reloss);
   fChain->SetBranchAddress("fqpi0pcflg", fqpi0pcflg, &b_fqpi0pcflg);
   fChain->SetBranchAddress("fqpi0mom1", fqpi0mom1, &b_fqpi0mom1);
   fChain->SetBranchAddress("fqpi0mom2", fqpi0mom2, &b_fqpi0mom2);
   fChain->SetBranchAddress("fqpi0momtot", fqpi0momtot, &b_fqpi0momtot);
   fChain->SetBranchAddress("fqpi0dconv1", fqpi0dconv1, &b_fqpi0dconv1);
   fChain->SetBranchAddress("fqpi0dconv2", fqpi0dconv2, &b_fqpi0dconv2);
   fChain->SetBranchAddress("fqpi0t0", fqpi0t0, &b_fqpi0t0);
   fChain->SetBranchAddress("fqpi0totmu", fqpi0totmu, &b_fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", fqpi0nll, &b_fqpi0nll);
   fChain->SetBranchAddress("fqpi0mass", fqpi0mass, &b_fqpi0mass);
   fChain->SetBranchAddress("fqpi0photangle", fqpi0photangle, &b_fqpi0photangle);
   fChain->SetBranchAddress("fqpi0pos", fqpi0pos, &b_fqpi0pos);
   fChain->SetBranchAddress("fqpi0dir1", fqpi0dir1, &b_fqpi0dir1);
   fChain->SetBranchAddress("fqpi0dir2", fqpi0dir2, &b_fqpi0dir2);
   fChain->SetBranchAddress("fqpi0dirtot", fqpi0dirtot, &b_fqpi0dirtot);
   fChain->SetBranchAddress("fqnmrfit", &fqnmrfit, &b_fqnmrfit);
   fChain->SetBranchAddress("fqmrifit", fqmrifit, &b_fqmrifit);
   fChain->SetBranchAddress("fqmrnring", fqmrnring, &b_fqmrnring);
   fChain->SetBranchAddress("fqmrpcflg", fqmrpcflg, &b_fqmrpcflg);
   fChain->SetBranchAddress("fqmrnll", fqmrnll, &b_fqmrnll);
   fChain->SetBranchAddress("fqmrtotmu", fqmrtotmu, &b_fqmrtotmu);
   fChain->SetBranchAddress("fqmrpid", fqmrpid, &b_fqmrpid);
   fChain->SetBranchAddress("fqmrmom", fqmrmom, &b_fqmrmom);
   fChain->SetBranchAddress("fqmrdconv", fqmrdconv, &b_fqmrdconv);
   fChain->SetBranchAddress("fqmreloss", fqmreloss, &b_fqmreloss);
   fChain->SetBranchAddress("fqmrt0", fqmrt0, &b_fqmrt0);
   fChain->SetBranchAddress("fqmrpos", fqmrpos, &b_fqmrpos);
   fChain->SetBranchAddress("fqmrdir", fqmrdir, &b_fqmrdir);
   fChain->SetBranchAddress("fqmsnfit", &fqmsnfit, &b_fqmsnfit);
   fChain->SetBranchAddress("fqmspcflg", fqmspcflg, &b_fqmspcflg);
   fChain->SetBranchAddress("fqmsnseg", fqmsnseg, &b_fqmsnseg);
   fChain->SetBranchAddress("fqmspid", fqmspid, &b_fqmspid);
   fChain->SetBranchAddress("fqmsifit", fqmsifit, &b_fqmsifit);
   fChain->SetBranchAddress("fqmsimer", fqmsimer, &b_fqmsimer);
   fChain->SetBranchAddress("fqmstotmu", fqmstotmu, &b_fqmstotmu);
   fChain->SetBranchAddress("fqmsnll", fqmsnll, &b_fqmsnll);
   fChain->SetBranchAddress("fqmsmom", fqmsmom, &b_fqmsmom);
   fChain->SetBranchAddress("fqmseloss", fqmseloss, &b_fqmseloss);
   fChain->SetBranchAddress("fqmst0", fqmst0, &b_fqmst0);
   fChain->SetBranchAddress("fqmspos", fqmspos, &b_fqmspos);
   fChain->SetBranchAddress("fqmsdir", fqmsdir, &b_fqmsdir);
   fChain->SetBranchAddress("prid_nused", &prid_nused, &b_prid_nused);
   fChain->SetBranchAddress("nused2", &nused2, &b_nused2);
   fChain->SetBranchAddress("t2knuefs", &t2knuefs, &b_t2knuefs);
   fChain->SetBranchAddress("t2knuenpi0", &t2knuenpi0, &b_t2knuenpi0);
   fChain->SetBranchAddress("t2knueflag", &t2knueflag, &b_t2knueflag);
   fChain->SetBranchAddress("totwgt", &totwgt, &b_totwgt);
   }

   Notify();
}

Bool_t fqEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fqEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fqEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef fqEvent_cxx