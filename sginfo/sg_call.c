#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#include "sginfo.h"
#include "sg_call.h"


static int str_ibegin(const char *s1, const char *s2)    // check whether s2 include s1
{
  char     u1, u2;

  while (*s1 && *s2)
  {
    u1 = toupper(*s1++);
    u2 = toupper(*s2++);
    if      (u1 < u2) return -1;
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;
  return 0;
}


int BuildSgInfo(T_SgInfo *SgInfo, const char *SgName)    //  build Sginfo
{
  int                VolLetter;
  const T_TabSgName  *tsgn;

  while (*SgName && isspace(*SgName)) SgName++;          // loop *SgName, skip blank, until find str of SgName

  VolLetter = 'A';

  tsgn = NULL;                                           // init the address of tsgn to NULL

  if (VolLetter)
  {
    tsgn = FindTabSgNameEntry(SgName, VolLetter);        // Run FindTabSgNameEntry() func to locate address of tsgn
    if (tsgn == NULL) return -1;                         /* no matching table entry */
    SgName = tsgn->HallSymbol;                           // Find Hallsymbol of the space group
  }

  /* Allocate memory for the list of Seitz matrices and
     a supporting list which holds the characteristics of
     the rotation parts of the Seitz matrices
   */

  SgInfo->MaxList = 192;                                 /* absolute maximum number of symops */

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));        // apply memory for ListSeitzMx

  if (SgInfo->ListSeitzMx == NULL) {                     // The address NULL　means is failed to apply memory
    SetSgError("Not enough core");
    return -1;
  }

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));      // apply memory for ListRotMxInfo

  if (SgInfo->ListRotMxInfo == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  InitSgInfo(SgInfo);
  SgInfo->TabSgName = tsgn;                               // assign tsgn to Sginfo, since we obtained it

  /* Translate the Hall symbol and generate the whole group
   */

  ParseHallSymbol(SgName, SgInfo);
  if (SgError != NULL) return -1;

  /* Do some book-keeping and derive crystal system, point group,
     and - if not already set - find the entry in the internal
     table of space group symbols
   */
  return CompleteSgInfo(SgInfo);                         // Finally, Run the CompleteSgInfo() func and return
}


T_SgInfo get_sginfo(char *spg)                         // function to get sginfo struct
{
  T_SgInfo  SgInfo;                                      // create instance of SgInfo struct

  if (BuildSgInfo(&SgInfo, spg) != 0)                    // build SgInfo
    fprintf(stderr, "%s\n", SgError);
  else
  {

    if (SgError)                                         // if SgError is not None, print SgError to stderr
      fprintf(stderr, "%s\n", SgError);
    }

  return SgInfo;
}



const T_TabSgName *get_tsgn(char *spg)                         // function to get sginfo struct
{
  int                VolLetter;
  const T_TabSgName  *tsgn;

  VolLetter = 'A';
  tsgn = NULL;                                          // init the address of tsgn to NULL

  if (VolLetter)
  {
    tsgn = FindTabSgNameEntry(spg, VolLetter);          // Run FindTabSgNameEntry() func to locate address of tsgn
  }

  return tsgn;
}



const char *get_laue(char *spg)                         // function to get laue group from space group
{
  T_SgInfo  SgInfo;                                      // create instance of SgInfo struct
  int pg_index;
  const char *laue_group;

  if (BuildSgInfo(&SgInfo, spg) != 0)                    // build SgInfo
    fprintf(stderr, "%s\n", SgError);
  else
  {
    pg_index = PG_Index(SgInfo.PointGroup);
    laue_group = PG_Names[PG_Index(LG_Code_of_PG_Index[pg_index])];

    if (SgError)                                         // if SgError is not None, print SgError to stderr
      fprintf(stderr, "%s\n", SgError);
    }

  return laue_group;
}

const char *get_pg(char *spg)                            // function to get point group from space group
{
  T_SgInfo  SgInfo;                                      // create instance of SgInfo struct
  int pg_index;
  const char *point_group;

  if (BuildSgInfo(&SgInfo, spg) != 0)                    // build SgInfo
    fprintf(stderr, "%s\n", SgError);
  else
  {

    pg_index = PG_Index(SgInfo.PointGroup);
    point_group = PG_Names[pg_index];

    if (SgError)                                         // if SgError is not None, print SgError to stderr
      fprintf(stderr, "%s\n", SgError);
    }

  return point_group;
}


const char *get_crystal_system(char *spg)
{
  T_SgInfo  SgInfo;
  const char *crystal_system;

  if (BuildSgInfo(&SgInfo, spg) != 0)
    fprintf(stderr, "%s\n", SgError);
  else
  {

    crystal_system = XS_Name[SgInfo.XtalSystem];

    if (SgError)
      fprintf(stderr, "%s\n", SgError);
    }

  return crystal_system;
}


int get_unique_axis(char *spg)
{
  T_SgInfo  SgInfo;
  int unique_axis;

  if (BuildSgInfo(&SgInfo, spg) != 0)
    fprintf(stderr, "%s\n", SgError);
  else
  {

    unique_axis = SgInfo.UniqueRefAxis;

    if (SgError)
      fprintf(stderr, "%s\n", SgError);
    }

  return unique_axis;
}


/* function to generate unique hkl list */
int hkl_list(char *spg, char *filename, int Maxh, int Maxk, int Maxl)
{
  int ListSysAbsent = 1;
  T_SgInfo SgInfo;

  if (BuildSgInfo(&SgInfo, spg) != 0)
    fprintf(stderr, "%s\n", SgError);

  int        h, k, l, iList, restriction, M, n, i;
  int        Minh, Mink, Minl;
  int        uvw[3];
  int        CCMx_PL[9], deterCCMx_LP = 0, hP, kP, lP;

  if (SgInfo.LatticeInfo->Code != 'P')
  {
    deterCCMx_LP = deterRotMx(SgInfo.CCMx_LP);
                 InverseRotMx(SgInfo.CCMx_LP, CCMx_PL);

    if (deterCCMx_LP < 1)
      return -1;
  }

  SetListMin_hkl(&SgInfo, Maxk, Maxl, &Minh, &Mink, &Minl);

  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Failed to open file");
    return 0;
  }

  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(&SgInfo, h, k, l, &restriction);
    if (SgError != NULL)
      return -1;

    M = BuildEq_hkl(&SgInfo, NULL, h, k, l);
    if (SgError != NULL)
      return -1;
    if (h==0 && k==0 && l==0) continue;

    if (iList == 0)
    {
      if ((iList = IsSuppressed_hkl(&SgInfo, Minh, Mink, Minl,
                                                  Maxk, Maxl,
                                               h,    k,    l)) != 0)
        continue;
      else
        fprintf(fp, "  %3d %3d %3d\n",
                            h, k, l);
    }
    else if (ListSysAbsent)
      continue;

}
  fclose(fp);
  return 0;
}


T_RTMx xyz2RTMx(char *SymXYZ)
{
  T_RTMx  SeitzMx;
  int        FacTr;

  FacTr = STBF;

  ParseSymXYZ(SymXYZ, &SeitzMx, FacTr);


  return SeitzMx;
}

// function to check whether  a hkl indice is systematically absent
int check_hkl_absent(char *spg, int h, int k, int l)              // TODO: revise me later !!!
{
  int *TH_Restriction;
  int absent;
  T_SgInfo SgInfo;
  TH_Restriction = NULL;

  if (BuildSgInfo(&SgInfo, spg) != 0)
    fprintf(stderr, "%s\n", SgError);

  else{
    absent = IsSysAbsent_hkl(&SgInfo, h, k, l, TH_Restriction);
  }

  return absent;
}

// function to get equivalent hkl strcut
T_Eq_hkl get_Eq_hkl(char *spg, int h, int k, int l)              // TODO: revise me later !!!
{
  T_Eq_hkl eq_hkl;
  T_SgInfo SgInfo;

  if (BuildSgInfo(&SgInfo, spg) != 0)
    fprintf(stderr, "%s\n", SgError);

  else{
    BuildEq_hkl(&SgInfo, &eq_hkl, h, k, l);
  }

  return eq_hkl;
}
