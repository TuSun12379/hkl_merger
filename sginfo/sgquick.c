#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#include "sginfo.h"

// function to check whether *s2 are contained in the first part of *s1
static int str_ibegin(const char *s1, const char *s2)
                                                            // function defined with static can only be used in this file. const keyword defined a paramter, the value of which can not be changed
{
  char     u1, u2;

  while (*s1 && *s2)                                        // loop *s1 and *s2, until one of them is empty
  {
    u1 = toupper(*s1++);                                    // toupper() translate lowercase letter into upper, first to translate, then to ++
    u2 = toupper(*s2++);
    if      (u1 < u2) return -1;                            // when there is different chr, just return -1 or 1
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;                                       // if *s2 "VolA" was not finished to loop (*s1 finished at the same time), which means the char in the first part of *s1 do not contain *s2
  return 0;                                                 // if nothing wrong,just return 0
}


int BuildSgInfo(T_SgInfo *SgInfo, const char *SgName)       // both the parameters are pointer, struct pointer or char pointer
{
  int                VolLetter;                             // used to represent convention of international table
  const T_TabSgName  *tsgn;                                 // T_TabSgName struct pointer, used to store information of space group after Hallsymbol was found


  /* look for "VolA", "VolI", or "Hall"
   */
// remove the blank in the *SgName until first chr was found
  while (*SgName && isspace(*SgName)) SgName++;             // check whether there is str in *Sgname，and isspace() is used to check whether there is blank，return 1 if yes.

  VolLetter = -1;                                           // init the volLetter as -1, used to define the kind of space group

  if      (isdigit(*SgName))                                // check whether the first chr is digit, if yes assign VolLetter as VolA
    VolLetter = 'A';     // 'A' will be translate into int, 65
  else if (str_ibegin(SgName, "VolA") == 0)                 // check whether "VolA" is in the first part of *SgName, if yes, return VolLetter as 'A'
  {
    VolLetter = 'A';
    SgName += 4;                                            // move the pointer to next address, but why to do this?
  }
  else if (   str_ibegin(SgName, "VolI") == 0
           || str_ibegin(SgName, "Vol1") == 0)              // check wheter "VolI" or "Vol1" are contained in the first part of *SgName
  {
    VolLetter = 'I';
    SgName += 4;
  }
  else if (str_ibegin(SgName, "Hall") == 0)                 // check wheter "Hall" is contained in the first part of *SgName
  {
    VolLetter = 0;                                          // return 0, if HallSymbol is found
    SgName += 4;
  }

  while (*SgName && isspace(*SgName)) SgName++;             // remove blank，until first char was found

  /* default is "VolA"
   */

  if (VolLetter == -1)                                      // use the default value, if VolLetter hasn't been changed.
    VolLetter = 'A';

  /* if we don't have a Hall symbol do a table look-up
   */

  tsgn = NULL;                                              // assign NULL to tsgn

/*if VolLetter exist, try to use FindTabSgNameEntry() to find the space group，
the parameters are SgName (pointer has been moved to next chr) VolLetter which represent the kind of space goup*/
  if (VolLetter)
  {
    tsgn = FindTabSgNameEntry(SgName, VolLetter);          // find the tsgn struct with FindTabSgNameEntry() function
    if (tsgn == NULL) return -1;                           /* no matching table entry */       // return -1, if nothing was found
    SgName = tsgn->HallSymbol;                             // if T_TabSgName was successfully constructed, HallSymbol was already stored in it
  }

  /* Allocate memory for the list of Seitz matrices and
     a supporting list which holds the characteristics of
     the rotation parts of the Seitz matrices
   */
// in theory, the maxmum number of Seitz matrices is 192, which can represent all point symmetry operations
  SgInfo->MaxList = 192;                                   /* absolute maximum number of symops */  // assign value to MaxList element of T_SgInfo struct

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));        // assign malloc for LsitSeitzMx in SgInfo struct

  if (SgInfo->ListSeitzMx == NULL) {                                  // if NULL return, means memory is not enough, just return -1
    SetSgError("Not enough core");
    return -1;
  }

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));      // assign malloc for ListRotMxInfo in SgInfo struct

  if (SgInfo->ListRotMxInfo == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  /* Initialize the SgInfo structure
   */

  InitSgInfo(SgInfo);                                                 // Init Sginfo struct, get all unknown elements a init value
  SgInfo->TabSgName = tsgn; /* in case we know the table entry */     // assign the tsgn (address of HallSymbol) to TabSgName element of SgInfo

  /* Translate the Hall symbol and generate the whole group
   */

  ParseHallSymbol(SgName, SgInfo);                                    // parse the HallSymbol, get all symmetry operation matrix
  if (SgError != NULL) return -1;                                     // if the standard error output is not NULL, then return -1.

  /* Do some book-keeping and derive crystal system, point group,
     and - if not already set - find the entry in the internal
     table of space group symbols
   */

  return CompleteSgInfo(SgInfo);                                      // Complete SgInfo
}


int main(int argc, char *argv[])                            // command line parameters, argc record the number of parameters
{
  T_SgInfo  SgInfo;                                         // create a T_SgInfo struct
  int iList;                                                // iList used to record the index which related to point group


  if (argc == 2)                                            // when there are two parameters
  {
    if (BuildSgInfo(&SgInfo, argv[1]) != 0)                 // parse the command line arg，search related space group，if not normally return, print error
      fprintf(stderr, "%s\n", SgError);
    else
    {
      ListSgInfo(&SgInfo, 1, 0, stdout);                    // print the result of ListSgInfo, if no error

      // iList = PG_Index(SgInfo.PointGroup);               // it looks that iList is a pointer
      printf("The index is: %n\n", iList);
      printf("The related point group is: %s\n",PG_Names[iList]);                                        // print point group of the defined space group
      printf("The related point group number is: %d\n",PG_Number(SgInfo.PointGroup));                    // print the number point group
      printf("The related laue group is:%s\n",PG_Names[PG_Index(LG_Code_of_PG_Index[iList])]);           // print related laue group, successful! I know that's very strange
      printf("The related space group is:%s\n",EI_Name[SgInfo.ExtraInfo]);                               // print related space group, failed !

      if (SgError)                                           //  if any error happend, return the error
        fprintf(stderr, "%s\n", SgError);
    }
  }

  return 0;  // 正常输出
}
