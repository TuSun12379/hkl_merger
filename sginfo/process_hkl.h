#include <stdio.h>
#define PI 3.14159265

typedef struct
  {
    int  *h;
    int  *k;
    int  *l;
    double *inte;
    double *sigma;
    int len;
  }
  RefList;

typedef struct {
  	// int have_parameters;
  	/* Crystallographic representation */
  	double a;	/* m */
  	double b;	/* m */
  	double c;	/* m */
  	double alpha;	/* Radians */
  	double beta;	/* Radians */
  	double gamma;	/* Radians */
  }
  UnitCell;

enum data_type
{
  EDT,
  XDS_SHELX,
  CAP,
  PETS2
};

enum data_type_out
{
  SHELX,
  GENERAL,
  SCALED
};

enum Laue_Asym_type
{
  ASYM1,
  ASYM2B,
  ASYM3,
  ASYM4A,
  ASYM4B,
  ASYM5C,
  ASYM5D,
  ASYM5E,
  ASYM6A,
  ASYM6B,
  ASYM7A,
  ASYM7B
};

int read_hkl(const char *filename, RefList *refl, enum data_type data_type);

void initRefList(RefList *refl);

int cell_get_parameters(const UnitCell *cell, double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma);

double resolution(UnitCell *cell, int h, int k, int l);

int write_hkl(RefList *refl, char *fname, enum data_type_out out);

int hkl2unique(const char *fn,  char *fn_out, char *spg);

int asym_1(int h, int k, int l);

int asym_2a(int h, int k, int l);

int asym_2b(int h, int k, int l);

int asym_2c(int h, int k, int l);

int asym_3(int h, int k, int l);

int asym_4a(int h, int k, int l);

int asym_4b(int h, int k, int l);

int asym_5a(int h, int k, int l);

int asym_5b(int h, int k, int l);

int asym_5c(int h, int k, int l);

int asym_5d(int h, int k, int l);

int asym_5e(int h, int k, int l);

int asym_6a(int h, int k, int l);

int asym_6b(int h, int k, int l);

int asym_7a(int h, int k, int l);

int asym_7b(int h, int k, int l);

enum Laue_Asym_type laue_type(char *spg);

int scale_two_file_I1overI2(const char *fn1, const char *fn2, char *fn_scale, UnitCell cell);

int I1overI2_scale(RefList *list1, RefList *list2, UnitCell *cell);
