#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "process_hkl.h"
#include "sg_call.h"
#include "sginfo.h"

int pow_int(int a, int b) {
  int i = 0;
  int res = 1;

  if (b < 0) printf("Wrong value of b!\n");

  if (b == 0) res = 1;

  else {
    for (i=0; i<b; i++) {
    res *= a;
  }
  }

  return res;
}

int asym_1(int h, int k, int l) {
  int ret = 0;

  if (h >= 0) return 1;
  if ((k < 0) && (h == 0)) ret = 0;
  if ((l < 0) && (h == 0) && (k == 0)) ret = 0;

  return ret;
}

int asym_2a(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (l >= 0)) ret = 1;
  if ((l < 0) && (k == 0)) ret = 0;

  return ret;
}

int asym_2b(int h, int k, int l) {
  int ret = 0;

  if ((k >= 0) && (l >= 0)) ret = 1;
  if ((h < 0) && (l == 0)) ret = 0;

  return ret;
}

int asym_2c(int h, int k, int l) {
  int ret = 0;

  if ((k >= 0) && (l >= 0)) ret = 1;
  if ((h < 0) && (k == 0)) ret = 0;

  return ret;
}

int asym_3(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= 0) && (l >= 0)) ret = 1;

  return ret;
}

int asym_4a(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= 0) && (l >= 0)) ret = 1;
  if ((k > 0) && (h == 0)) ret = 0;

  return ret;
}

int asym_4b(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= h) && (l >= 0)) ret = 1;

  return ret;
}

int asym_5a(int h, int k, int l) {
  int ret = 0;

  if ((k >= 0) && (l >= 0) && (k >= h) && (l >= h)) ret = 1;
  if ((k > h) && (h == l)) ret = 0;
  if ((l < 0) && (h == 0)) ret = 0;

  return ret;
}

int asym_5b(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= 0))  ret = 1;
  if ((l < 0) && (k == 0)) ret = 0;
  if ((l <= 0) && (h == 0)) ret = 0;

  return ret;
}

int asym_5c(int h, int k, int l) {
  int ret = 0;

  if ((k >= 0) && (l >= k) && (k >= h))  ret = 1;
  if ((h < 0) && (abs(h) > l) && (l == 0)) ret = 1;

  return ret;
}

int asym_5d(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= h))  ret = 1;
  if ((l < 0) && (h == 0)) ret = 0;

  return ret;
}

int asym_5e(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= 0) && (l >= 0))  ret = 1;
  if ((h < k) && (l == 0)) ret = 0;

  return ret;
}

int asym_6a(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= 0) && (l >= 0))  ret = 1;
  if ((k > 0) && (h == 0)) ret = 0;

  return ret;
}

int asym_6b(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= h) && (l >= 0))  ret = 1;

  return ret;
}

int asym_7a(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= h) && (l >= h))  ret = 1;
  if ((k > h) && (h == l)) ret = 0;

  return ret;
}

int asym_7b(int h, int k, int l) {
  int ret = 0;

  if ((h >= 0) && (k >= h) && (l >= k))  ret = 1;

  return ret;
}

enum Laue_Asym_type laue_type(char *spg) {

  const char *laue;
  char *asym1 = "-1";
  char *asym2b = "2/m";
  char *asym3 = "mmm";
  char *asym4a = "4/m";
  char *asym4b = "4/mmm";
  char *asym5c = "-3";
  char *asym5d = "-31m";
  char *asym5e = "-3m1";
  char *asym6a = "6/m";
  char *asym6b = "6/mmm";
  char *asym7a = "m-3";
  char *asym7b = "m-3m";
  enum Laue_Asym_type asym_type;

  laue = get_laue(spg);
  printf("The laue group of %s is %s.\n", spg, laue);

  if (strcmp(laue, asym1) == 0) asym_type = ASYM1;
  else if (strcmp(laue, asym2b) == 0) asym_type = ASYM2B;
  else if (strcmp(laue, asym3) == 0) asym_type = ASYM3;
  else if (strcmp(laue, asym4a) == 0) asym_type = ASYM4A;
  else if (strcmp(laue, asym4b) == 0) asym_type = ASYM4B;
  else if (strcmp(laue, asym5c) == 0) asym_type = ASYM5C;
  else if (strcmp(laue, asym5d) == 0) asym_type = ASYM5D;
  else if (strcmp(laue, asym5e) == 0) asym_type = ASYM5E;
  else if (strcmp(laue, asym6a) == 0) asym_type = ASYM6A;
  else if (strcmp(laue, asym6b) == 0) asym_type = ASYM6B;
  else if (strcmp(laue, asym7a) == 0) asym_type = ASYM7A;
  else if (strcmp(laue, asym7b) == 0) asym_type = ASYM7B;
  else printf("Laue group determine error!.\n");

  return asym_type;

}

int hkl2unique(const char *fn, char *fn_out, char *spg) {

  RefList refl;
  int len;
  int i, j, z;
  enum Laue_Asym_type asym_type;
  enum data_type data_type = XDS_SHELX;
  enum data_type_out out = SCALED;

  asym_type = laue_type(spg);

  if (read_hkl(fn, &refl, data_type) != 0) printf("Reflection list resording error!.\n");

  len = refl.len;

  for(i=0; i<len; i++) {
    T_Eq_hkl Eq_hkl;
    int N;
    int h, k, l;
    int rest;

    Eq_hkl = get_Eq_hkl(spg, refl.h[i], refl.k[i],refl.l[i]);
    N = Eq_hkl.N;

    for (j=0; j<N; j++) {

      h = Eq_hkl.h[j];
      k = Eq_hkl.k[j];
      l = Eq_hkl.l[j];

    for (z=0; z<2; z++) {
      h =  h * pow_int(-1, z);
      k =  k * pow_int(-1, z);
      l =  l * pow_int(-1, z);

      switch (asym_type) {
      case ASYM1:
        rest = asym_1(h, k, l);
        break;

      case ASYM2B:
        rest = asym_2b(h, k, l);
        break;

      case ASYM3:
        rest = asym_3(h, k, l);
        break;

      case ASYM4A:
        rest = asym_4a(h, k, l);
        break;

      case ASYM4B:
        rest = asym_4b(h, k, l);
        break;

      case ASYM5C:
        rest = asym_5c(h, k, l);
        break;

      case ASYM5D:
        rest = asym_5d(h, k, l);
        break;

      case ASYM5E:
        rest = asym_5e(h, k, l);
        break;

      case ASYM6A:
        rest = asym_6a(h, k, l);
        break;

      case ASYM6B:
        rest = asym_6b(h, k, l);
        break;

      case ASYM7A:
        rest = asym_7a(h, k, l);
        break;

      case ASYM7B:
        rest = asym_7b(h, k, l);
        break;
      }

      if (rest == 1) {
        refl.h[i] = h;
        refl.k[i] = k;
        refl.l[i] = l;
        // printf("%4d%4d%4d\n", h, k, l);
      }

    }


    }
    // printf("\n");

    // printf("The equivalent reflecion number of %4d %4d %4d is %d.\n", refl.h[i], refl.k[i],refl.l[i], 2*N);
  }

  write_hkl(&refl, fn_out, out);

  return 0;
}

void initRefList(RefList *refl) {
    refl->h = NULL;
    refl->k = NULL;
    refl->l = NULL;
    refl->inte = NULL;
    refl->sigma = NULL;
    refl->len = 0;
}

int read_hkl(const char *filename, RefList *refl, enum data_type data_type)
{
  FILE *fp;
  int max_n=1024;

  initRefList(refl);

	if ( filename == NULL ) {
		printf("Can't open file %s.\n", filename);
    return 1;
	} else {
		fp = fopen(filename, "r");
	}

  if(fp==NULL) {
    printf("Can't open file %s.\n", filename);
    return 1; }
    else {
    printf("The file: %s is successfully read!\n", filename);
  }

  char line[1024];
  char *rval = NULL;
  int i = 0;
  rval = fgets(line, 1023, fp);
	if ( rval == NULL ) {
    printf("Can not read the first line !\n");
    return 1;
  }

if (refl==NULL) printf("refl is NULL!");
  refl->h = malloc(max_n*sizeof(int));
  refl->k = malloc(max_n*sizeof(int));
  refl->l = malloc(max_n*sizeof(int));
  refl->inte = malloc(max_n*sizeof(double));
  refl->sigma = malloc(max_n*sizeof(double));

  if ((refl->h==NULL) || (refl->k==NULL) || (refl->l==NULL) ||
            (refl->inte==NULL) || (refl->sigma==NULL)) {
		printf("Failed to allocate memory for recording data into reflection list.\n");
		return 1; }

  while ( rval != NULL )
  {
    int h, k, l;
    float intensity, sigma;
    int r;

    switch (data_type) {
      case EDT :
      r = sscanf(line, "%i %i %i %f %f", &h, &k, &l, &intensity, &sigma);
      break;

      case XDS_SHELX :
      r = sscanf(line, "%i %i %i %f %f", &h, &k, &l, &intensity, &sigma);
      break;

      case CAP :
      r = sscanf(line, "%i %i %i %f %f", &h, &k, &l, &intensity, &sigma);
      break;              // revise me later !!!

      case PETS2 :
      r = sscanf(line, "%i %i %i %f %f", &h, &k, &l, &intensity, &sigma);
      break;
    }

    if ( r != 5 ) {
      printf("Bad line '%s'\n", line);
    }

    if (h==0 && k==0 && l==0) {
      printf("Bad line: %4d%4d%4d; it will not be recorded.\n", h, k, l);
    } else{

      refl->h[i] = h;
      refl->k[i] = k;
      refl->l[i] = l;
      refl->inte[i] = intensity;
      refl->sigma[i] = sigma;
      i++;
    }

    if ( i == max_n ){
      max_n += 1024;
      refl->h = realloc(refl->h, max_n*sizeof(int));
      refl->k = realloc(refl->k, max_n*sizeof(int));
      refl->l = realloc(refl->l, max_n*sizeof(int));
      refl->inte = realloc(refl->inte, max_n*sizeof(double));
      refl->sigma = realloc(refl->sigma, max_n*sizeof(double));

    if ( (refl->h==NULL) || (refl->k==NULL) || (refl->l==NULL) ||
      (refl->inte==NULL) || (refl->sigma==NULL)) {
      printf("Failed to allocate more memory for recording data.\n");
      return 1; }
      }

    rval = fgets(line, 1023, fp);
  }

  refl->len = i;
  fclose(fp);

return 0;
}

int cell_get_parameters(const UnitCell *cell, double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma)
{
		/* Direct response */
		*a = cell->a;
		*b = cell->b;
		*c = cell->c;
		*alpha = cell->alpha;
		*beta = cell->beta;
		*gamma = cell->gamma;

	return 0;
}

/* Return sin(theta)/lambda = 1/2d.  Multiply by two if you want 1/d */
double resolution(UnitCell *cell, int h, int k, int l)
{
	double a, b, c, alpha, beta, gamma;

	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

	const double Vsq = a*a*b*b*c*c*(1 - cos(alpha)*cos(alpha)
	                                  - cos(beta)*cos(beta)
	                                  - cos(gamma)*cos(gamma)
				          + 2*cos(alpha)*cos(beta)*cos(gamma) );

	const double S11 = b*b*c*c*sin(alpha)*sin(alpha);
	const double S22 = a*a*c*c*sin(beta)*sin(beta);
	const double S33 = a*a*b*b*sin(gamma)*sin(gamma);
	const double S12 = a*b*c*c*(cos(alpha)*cos(beta) - cos(gamma));
	const double S23 = a*a*b*c*(cos(beta)*cos(gamma) - cos(alpha));
	const double S13 = a*b*b*c*(cos(gamma)*cos(alpha) - cos(beta));

	const double brackets = S11*h*h + S22*k*k + S33*l*l
	                         + 2*S12*h*k + 2*S23*k*l + 2*S13*h*l;
	const double oneoverdsq = brackets / Vsq;
	const double oneoverd = sqrt(oneoverdsq);

	return oneoverd / 2;
}

double d_value(UnitCell *cell, int h, int k, int l)
{ double d;

  d = resolution(cell, h, k, l);
  d = 0.5*(1/d);

	return d;
}

int I1overI2_scale(RefList *list1, RefList *list2, UnitCell *cell)
{
  int len1, len2;
  len1 = list1->len;
  len2 = list2->len;

  printf("I1/I2 scaling process started.\n");
  printf("The length of list1: %d and list2: %d.\n", len1, len2);

  double s;
	double top = 0.0;
	double bot = 0.0;
  int i, j;
  int n = 0;


  for(i=0; i<len1; i++)
	{
		double I1, I2;
		double res;
		for(j=0; j<len2; j++)
		{
			if (list1->h[i]==list2->h[j] && list1->k[i]==list2->k[j] && list1->l[i]==list2->l[j])
			{

				I1 = list1->inte[i];
			  I2 = list2->inte[j];

          if ( (I1 <= 0.0) || (I2 <= 0.0)) {
            printf("Negative value will be ignored.\n");
            continue;
            }

          else {
			    top += I1;
          bot += I2;
          n++;
          }
			}
		}
  }

s = top / bot;

printf("The same reflection number is: %d.\n", n);

// check for number of reflections used for scaling
if ( n < 2 ) {
  printf("Not enough reflections for scaling\n");
  return 1;
}

// loop the list2 and scale it.
for(i=0; i<len2; i++)
{

  double inte_temp;
  double sigma_temp;

  inte_temp = list2->inte[i];
  sigma_temp = list2->sigma[i];



  list2->inte[i] = inte_temp*s;
  list2->sigma[i] = sigma_temp*s;
}

  printf("general wilson scale value: %8.4f.\n", s);

  return 0;
}

int scale_two_file_I1overI2(const char *fn1, const char *fn2, char *fn_scale, UnitCell cell) {
  RefList refl1;
  RefList refl2;
  enum data_type data_type = XDS_SHELX;
  enum data_type_out out = SCALED;

  if (read_hkl(fn1, &refl1, data_type) != 0) printf("Reflection list resording error!.\n");
  if (read_hkl(fn2, &refl2, data_type) != 0) printf("Reflection list resording error!.\n");

  I1overI2_scale(&refl1, &refl2, &cell);

  // write_hkl(&refl1, fn3, out);
  write_hkl(&refl2, fn_scale, out);

return 0;
}

int write_hkl(RefList *refl, char *fname, enum data_type_out out)
{
	int len;
  int i;
	FILE *fp;
  double max_inte = 0.0;

  len = refl->len;
  printf("%d reflections will be writen into file: %s.\n", len, fname);
	if((fp=fopen(fname,"w"))==NULL)
	{
		fprintf(stderr,"Can't open %s.\n",fname);
		// exit(1);
	}

  	for(i=0;i<len;i++)
	{
		if (refl->inte[i] > max_inte) {
      max_inte = refl->inte[i];
    }
	}
  // printf("The maximum intensity is: %f.\n", max_inte);

	for(i=0;i<len;i++)
	{
    double intensity;
    double sigma;

    if ( (refl->h[i] == 0) && (refl->k[i] == 0) && (refl->l[i] == 0) ) {
      printf("Zero indices 000 will be ignored.");
      continue;
    }

    switch (out)
    {
    case SHELX :
      intensity = refl->inte[i] * 99999 / max_inte;
      sigma = refl->sigma[i] * 99999 / max_inte;
      fprintf(fp, "%4d%4d%4d%8.2f%8.2f\n", refl->h[i], refl->k[i], refl->l[i], intensity, sigma);
      break;

    case GENERAL :
      intensity = refl->inte[i] * 99999 / max_inte;
      sigma = refl->sigma[i] * 99999 / max_inte;
      fprintf(fp, "%4d%4d%4d%12.2f%12.2f\n", refl->h[i], refl->k[i], refl->l[i], intensity, sigma);
      break;

    case SCALED :
      intensity = refl->inte[i];
      sigma = refl->sigma[i];
      // printf("Nomarlization will not be applied when write file: %s.\n",fname);
      fprintf(fp, "%4d%4d%4d%12.2f%12.2f\n", refl->h[i], refl->k[i], refl->l[i], intensity, sigma);
      break;
    }
	}

	fclose(fp);

return 0;
}
