#include <stdio.h>
#include "sginfo.h"

const char *get_laue(char *spg);                                           // call laue group from space group
int hkl_list(char *spg, char *filename, int Maxh, int Maxk, int Maxl);     // generate unique hkl list
const char *get_pg(char *spg);                                             // get pointer group
const char *get_crystal_system(char *spg);                                 // get crystal system
int get_unique_axis(char *spg);                                            // get unique axis
T_SgInfo get_sginfo(char *spg);                                            // get sginfo struct
const T_TabSgName *get_tsgn(char *spg);                                    // get tsgn struct
T_RTMx xyz2RTMx(char *SymXYZ);                                             // get RTMx with xyz symmetry
int check_hkl_absent(char *spg, int h, int k, int l);
T_Eq_hkl get_Eq_hkl(char *spg, int h, int k, int l);
