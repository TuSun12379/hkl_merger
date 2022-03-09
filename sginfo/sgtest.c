#include <stdio.h>
#include "sginfo.h"
#include "sg_call.h"

int main(void)
{
    char *spg  = "p6122";
    char *filename = "independent.txt";
    int Maxh = 50, Maxk = 50, Maxl = 50;
    hkl_list(spg, filename, Maxh, Maxk, Maxl);

    const char *laue_spg;
    laue_spg = get_laue(spg);
    printf("The laue group you wanted is: %s.\n", laue_spg);

    const char *point_group;
    point_group = get_pg(spg);
    printf("The point group you wanted is: %s.\n", point_group);

    const char *crystal_system;
    crystal_system = get_crystal_system(spg);
    printf("The crystal system you wanted is: %s.\n", crystal_system);

    int unique_axis;
    unique_axis = get_unique_axis(spg);
    printf("The uniuqe axis you wanted is: %c.\n", unique_axis);

    int pg_index;
    const char *laue_group;
    T_SgInfo  SgInfo;
    SgInfo = get_sginfo(spg);
    pg_index = PG_Index(SgInfo.PointGroup);
    laue_group = PG_Names[PG_Index(LG_Code_of_PG_Index[pg_index])];
    printf("The laue group you wanted is: %s.\n", laue_group);

    const T_TabSgName  *tsgn;
    int SgNumber;
    tsgn = get_tsgn(spg);
    SgNumber = tsgn->SgNumber;
    printf("The sg number you wanted is: %d.\n", SgNumber);

    T_RTMx rtms;
    char *SymXYZ = "-y, x-y, z+1/3";
    int a;
    rtms = xyz2RTMx(SymXYZ);
    a = rtms.a[0];
    printf("The first number of a is: %d.\n", a);

    int h=0, k=0, l=6;
    int absent;
    absent = check_hkl_absent(spg, h, k, l);
    printf("%d%d%d is absent?: %d.\n", h, k, l, absent);

    T_Eq_hkl eq_hkl;
    int M;
    eq_hkl = get_Eq_hkl(spg, h, k, l);
    M = eq_hkl.M;
    printf("The multiplicity of 112 is: %d.\n", M);

    return 0;
}
