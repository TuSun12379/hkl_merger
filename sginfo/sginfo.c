/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


/*
  Macintosh extras (Courtesy Jon Tischler <TischlerJZ@ornl.gov>)
 */
#if defined(__THINK__) || defined(__MWERKS__)
#include <console.h>
#define CONSOLE_LINES   36  /* number of lines to use for console */
#define CONSOLE_COLUMNS 90  /* number of columns to use for console */
#ifdef __MWERKS__
#include <sioux.h>
#endif
#endif


#define AppMalloc(ptr, n) (ptr) = malloc((n) * sizeof (*(ptr)))
#define AppFree(ptr, n) free(ptr)


#define SGCOREDEF__
#include "sginfo.h"


#if USE_GS_SI

static int PrimitiveRotMx(const int *CCMx_LP, int *RotMx, const int *CCMx_PL,
                          int deterCCMx_LP)
{
  int       i;
  int       BufMx[9];


  /* Mp = Tlp . Mz . Tpl */

  RotMxMultiply(BufMx, RotMx, CCMx_PL);
  RotMxMultiply(RotMx, CCMx_LP, BufMx);

  for (i = 0; i < 9; i++)
  {
    if (RotMx[i] % deterCCMx_LP) {
      SetSgError("Internal Error: PrimitiveRotMx()");
      return -1;
    }
  }

  for (i = 0; i < 9; i++)
    RotMx[i] /= deterCCMx_LP;

  return 0;
}


static int Find_si(T_SgInfo *SgInfo)
{
  static const int Tab_si_Vector[] =
    {
       1,  0,  0,   0, /*  h      */
       0,  1,  0,   1, /*  k      */
       0,  0,  1,   2, /*  l      */
       1,  1,  0,   0, /*  h+k    */
       1, -1,  0,   0, /*  h-k    */
       0,  1,  1,   1, /*  k+l    */
       0,  1, -1,   1, /*  k-l    */
       1,  0,  1,   1, /*  h+l    */
       1,  0, -1,   1, /*  h-l    */
       1,  1,  1,   0, /*  h+k+l  */
       1,  1, -1,   0, /*  h+k-l  */
       1, -1,  1,   0, /*  h-k+l  */
      -1,  1,  1,   0, /* -h+k+l  */
       2,  1, -1,   0, /*  2h+k-l */
       2, -1,  1,   0, /*  2h-k+l */
      -1,  2,  1,   0, /* -h+2k+l */
       1,  2, -1,   0, /*  h+2k-l */
      -1,  1,  2,   0, /* -h+k+2l */
       1, -1,  2,   0  /*  h-k+2l */
    };

  static int nTab_si_Vector
     = sizeof Tab_si_Vector / sizeof (*Tab_si_Vector) / 4;

  int        deterCCMx_LP, CCMx_PL[9];
  int        i, itabsiv;
  int        nLoopInv, iLoopInv, n_si_v, i_si_v;
  int        n, m, l;
  int        IsFine;
  int        item[3];
  int        R_I[9], si_Buf[9];
  int        iList;
  T_RTMx     *lsmx;
  const int  *tabsiv;


  if (SgInfo->LatticeInfo->Code != 'P')
  {
    deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP);
                 InverseRotMx(SgInfo->CCMx_LP, CCMx_PL);

    if (deterCCMx_LP < 1)
      goto ReturnError;
  }

  nLoopInv = Sg_nLoopInv(SgInfo);

  SgInfo->n_si_Vector = n_si_v = 0;

  for (i = 0; i < 9; i++)
    SgInfo->si_Vector[i] = 0;

  for (i = 0; i < 3; i++)
  {
    SgInfo->si_Modulus[i] = 1;
    item[i] = 1;
  }

  tabsiv = Tab_si_Vector;

  for (itabsiv = 0; itabsiv < nTab_si_Vector; itabsiv++, tabsiv += 4)
  {
    IsFine = 1;
    m = -1;

    for (iList = 0; IsFine && iList < SgInfo->nList; iList++)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      for (iLoopInv = 0; IsFine && iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv == 0)
          for (i = 0; i < 9; i++)
          {
            if (i % 4) R_I[i] =  lsmx->s.R[i];
            else       R_I[i] =  lsmx->s.R[i] - 1;
          }
        else
          for (i = 0; i < 9; i++)
          {
            if (i % 4) R_I[i] = -lsmx->s.R[i];
            else       R_I[i] = -lsmx->s.R[i] - 1;
          }

        if (SgInfo->LatticeInfo->Code != 'P')
        {
          if (PrimitiveRotMx(SgInfo->CCMx_LP, R_I, CCMx_PL,
                                deterCCMx_LP) < 0)
            return -1;
        }

        for (i = 0; IsFine && i < 3; i++)
        {
          n =  tabsiv[0] * R_I[i * 3 + 0];
          n += tabsiv[1] * R_I[i * 3 + 1];
          n += tabsiv[2] * R_I[i * 3 + 2];
          n = abs(n);

          if (n == 1)
            IsFine = 0;
          else if (m < 2)
            m = n;
          else if (n > 0 && n != m)
            IsFine = 0;
        }
      }
    }

    if (IsFine)
    {
#if DEBUG_Find_si
      fprintf(stdout, "H-Kt %2d %2d %2d   %d\n",
        tabsiv[0], tabsiv[1], tabsiv[2], m);
#endif

      l = tabsiv[3];

      while (item[l] > 1) /* just "if", see break's */
      {
        if (m == item[l]) break;

        if (m == 3 && (   SgInfo->XtalSystem != XS_Trigonal
                       || SgInfo->UniqueDirCode != '=')) break;

        if (m == 4 && (   SgInfo->XtalSystem == XS_Triclinic
                       || SgInfo->XtalSystem == XS_Monoclinic)) break;

        if (m == 2) break;

        /* if (m > 1 || m != 4) break; */

        n_si_v--;
        item[l] = 1;
        break;
      }

      if (item[l] == 1)
      {
        if (itabsiv > 12)
          n_si_v = 0;

        item[l] = m;
        SgInfo->si_Modulus[n_si_v] = m;

        n = n_si_v * 3;
        for (i = 0; i < 3; i++)
          SgInfo->si_Vector[n++] = tabsiv[i];

        n_si_v++;
      }
    }
  }

#if DEBUG_Find_si
  fprintf(stdout, "H-Kt\n");
#endif

  if (SgInfo->LatticeInfo->Code != 'P')
  {
#if DEBUG_Find_si
    for (i = 0; i < n_si_v; i++)
      fprintf(stdout, "H-Kp %2d %2d %2d   %d\n",
        SgInfo->si_Vector[i * 3 + 0],
        SgInfo->si_Vector[i * 3 + 1],
        SgInfo->si_Vector[i * 3 + 2],
        SgInfo->si_Modulus[i]);
    fprintf(stdout, "H-Kp\n");
#endif

    for (i_si_v = 0; i_si_v < n_si_v; i_si_v++)
    {
      for (i = 0; i < 3; i++)
      {
        si_Buf[i_si_v * 3 + i]
          =   SgInfo->si_Vector[i_si_v * 3 + 0] * CCMx_PL[i * 3 + 0]
            + SgInfo->si_Vector[i_si_v * 3 + 1] * CCMx_PL[i * 3 + 1]
            + SgInfo->si_Vector[i_si_v * 3 + 2] * CCMx_PL[i * 3 + 2];
      }
    }

    for (i = 0; i < i_si_v * 3; i++)
    {
      if (si_Buf[i] % deterCCMx_LP)
      {
        n = i / 3; n *= 3;
        fprintf(stdout, " %3d %3d %3d\n",
          si_Buf[n + 0], si_Buf[n + 1], si_Buf[n + 2]);
        goto ReturnError;
      }

      SgInfo->si_Vector[i] = si_Buf[i] / deterCCMx_LP;
    }
  }

  SgInfo->n_si_Vector = n_si_v;
  return n_si_v;

  ReturnError:

  SetSgError("Internal Error: Find_si()");
  return -1;
}


static int Try_GS_si(T_SgInfo *SgInfo)
{
  int       h, k, l, iList;
  int       Maxh, Maxk, Maxl;
  int       Minh, Mink, Minl;
  int       nTestField, *TestField;
  int       nProperty, *Property, *pp;
  int       IsFine, would_be, is;


  SgInfo->n_si_Vector = -1;

                       nTestField = 12 * 12 * 12;
  AppMalloc(TestField, nTestField);
  if (TestField == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  MarkLegalOrigins(SgInfo, TestField);

  Maxh = Maxk = Maxl = 7;
  SetListMin_hkl(SgInfo, Maxk, Maxl, &Minh, &Mink, &Minl);

  nProperty =   (Maxh - Minh + 1)
              * (Maxk - Mink + 1)
              * (Maxl - Minl + 1);

  AppMalloc(Property, nProperty);
  if (Property == NULL) {
    SetSgError("Not enough core");
    AppFree(TestField, nTestField);
    return -1;
  }

  pp = Property;
  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(SgInfo, h, k, l, NULL);
    if (SgError != NULL)
    {
      AppFree(Property, nProperty);
      AppFree(TestField, nTestField);
      return -1;
    }

    if (iList == 0)
      *pp++ = Verify_si(h, k, l, TestField);
    else
      *pp++ = -1;
  }

  if (Find_si(SgInfo) >= 0)
  {
    IsFine = 1;

    pp = Property;
    for (h = Minh; IsFine && h <= Maxh; h++)
    for (k = Mink; IsFine && k <= Maxk; k++)
    for (l = Minl; IsFine && l <= Maxl; l++)
    {
      is = *pp++;

      if (is >= 0)
      {
        would_be = Is_si(SgInfo, h, k, l);
        if (is != would_be)
          IsFine = 0;
      }
    }

    if (IsFine)
    {
      AppFree(Property, nProperty);
      AppFree(TestField, nTestField);
      return 0;
    }
  }

  SetSgError("Internal Error: Can't determine s.i. vectors and moduli");

  AppFree(Property, nProperty);
  AppFree(TestField, nTestField);

  return -1;
}

#endif /* USE_GS_SI */


static const char *progn = "sginfo";


static void progerror(const char *message)
{
  fflush(stdout);
  fprintf(stderr, "%s: %s\n", progn, message);
  exit(1);
}


static void NotEnoughCore(void)
{
  progerror("Not enough core");
}


static void PrintClearSgError(int ClearError, int CertainSgError)
{
  if (CertainSgError && SgError == NULL)
    SetSgError("Internal Error: SgError not set but should be");

  if (SgError)
  {
    fprintf(stdout, "%s: %s\n", progn, SgError);
    if (ClearError == 0) exit(1);
    SgError = NULL;
  }
}


static int str_icmp(const char *s, const char *t)  // used to check wheter two str are same or not; very similar to str_ibegin()
{
  char     cs, ct;

  while (*s || *t)  // 
  { cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }
  return 0;
}


static int str_ibegin(const char *s1, const char *s2) /* string ignore-case */
{                                                     /* begin              */
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

/*不可修改字符串数组, 只能在本文件中使用*/
static const char *LegendTabSgName[] =
  {
    "",
    "  Extensions",
    "  ----------",
    "    Monoclinic             unique axis b   unique axis c   unique axis a",
    "                             abc   c-ba      abc   ba-c      abc   -acb",
    "                            ------------    ------------    ------------",
    "             cell choice 1   :b1   :-b1      :c1   :-c1      :a1   :-a1",
    "                         2   :b2   :-b2      :c2   :-c2      :a2   :-a2",
    "                         3   :b3   :-b3      :c3   :-c3      :a3   :-a3",
    "",
    "    Orthorhombic   :ba-c    change of basis abc -> ba-c",
    "                   :1       origin choice 1",
    "                   :2ba-c   origin choice 2, change of basis abc -> ba-c",
    "",
    "    Tetragonal     :1       origin choice 1",
    "           Cubic   :2       origin choice 2",
    "",
    "    Trigonal       :H       hexagonal    axes",
    "                   :R       rhombohedral axes",
    "",
    "  Number   Schoenflies   Hermann-Mauguin             Hall",
    "  ------   -----------   ---------------             ----",
    NULL,
  };


static void ListTabSgName(int WantedSgNumber, int VolLetter, FILE *fpout)  
/*实参为F_ListTable, 以及VolLetter符号，第三个参数为标准输出流，是一个File类型指针*/
{
  int                i;
  const char         *sgl, *ext, **ltsgn;
  const T_TabSgName  *tsgn, *show, *show_later;  // tsgn, show和show_later都为T_TabSgName结构体指针


  if (WantedSgNumber == -1)  // 如果在最初的命令行解析中查找到了-ListTab字符串，但是其值被设定为了-1, 此时的设定是打印所有230个空间群的符号对照表
    for (ltsgn = LegendTabSgName; *ltsgn; ltsgn++)  // 说明没有找到对应空间群，遍历LegenTabSgName数组，只要对应字符串不为空，打印对应字符串元素
     fprintf(fpout, "%s\n", *ltsgn); // 首先将LegendTabSgName[]数组中的所有元素全部打印纸标准输出流

  if (VolLetter == '1') // 如果VolLetter为'1', 则将VolLetter设置为'I'。其它情况VolLetter应对应字母
    VolLetter = 'I';
  else
    VolLetter = toupper(VolLetter); // 当VolLetter不为1时，将字符串小写字符都转换为大写

  show = show_later = NULL;  // 暂时将show和show_later指针设为空值
/*TabSgName为T_TabSgName结构体数组指针，该数组中记录着所有230个空间群的各种字符对照表*/
  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++) // 遍历230个空间群表，查找和wantedSgNumber数字对应的空间群
  /*TabSgName[]是T_TabSgname结构体为模板声明的结构体数组，里面记录了从1到230号空间群对应的空间群字符以及VolA简写及完整VolA字符信息对照表。
  整个是一个结构体数组。将其首地址赋值给tsgn。遍历TabSgName[]数组，当对应的HallSymbol地址为空时跳出循环*/
  {
    if (   WantedSgNumber == -1
        || WantedSgNumber == tsgn->SgNumber) // F_Listtable为-1或当在空间群表中查找到对应的空间群时继续执行下面的操作
    {

/*#######################################  111111111111111111111111  #######################################################*/
      if (tsgn->SgNumber >= 3 && tsgn->SgNumber < 16) // 1.&&&&& 3到15号空间群为单斜
      {
        if (VolLetter == 'I')  /* ####1.当VolLetter为I时 ####*/
        {
               ext = tsgn->Extension;  // 在空间群表对应的结构体元素中查找完全空间群字符串
          if (*ext == '-')  // 如果第一个字符为'-'则地址前进一位
               ext++;

          if (       tsgn->Extension[0] == 'b'
              && (   tsgn->Extension[1] == '\0'
                  || tsgn->Extension[1] == '1')) // 当相应空间群字符串首字母为b，第二个字母为空或'1'时，将tsgn赋值给show_later
            show_later = tsgn; 
          else if (  ext[0] == 'c')  // 当第一个标记符号为 '-',空间群标记符号第二个字符为c时
          {
            if (ext[1] == '\0')  // 当第一个标记符号为 '-',空间群标记符号第二个字符为c, 第三个字母为空值时将tsgn赋值给show指针
              show = tsgn; // 将tsgn赋值给show
            else  // 当第一个标记符号为 '-',第二个字符为c,且第三个不为空值时，遍历 tsgn->SgLabels中的字符，如果存在两个'='，则将tsgn赋值给show指针
            {
              i = 0;
              for (sgl = tsgn->SgLabels; *sgl; sgl++)  
                if (*sgl == '=') i++; // 遍历 tsgn->SgLabels 中的所有字符，当*sgl等于'='时i自加1

              if (i == 2)  // 如果该字符串中存在两个'='，则将tsgn赋值给变量show
                show = tsgn;
            }
          }
        }
        else if (VolLetter == 'A') /* ####2.当VolLetter为A时 ####*/
        {
          if (   tsgn->Extension[0] != '-'
              && tsgn->Extension[0] != 'a')  // 当*tsgn结构体空间群标记符号首字母不为'-'或'a'时，将tsgn赋值给show
            show = tsgn;
        }
        else
          show = tsgn;  /* ####3.当VolLetter即不为I也不为A时，直接将tsgn赋值给show ####*/
      }

/*#######################################  2222222222222222222222  #######################################################*/
      else if (   tsgn->Extension[0] == 'H'  // 单斜晶系之外的情况,如果*tsgn结构体中空间群符号首字母为'H'同时VolLetteR为'I'时直接将tsgn赋值给show_later
               && VolLetter == 'I')
        show_later = tsgn;

/*#######################################  333333333333333333333  #######################################################*/
      else if (   VolLetter == 'A'   // 单斜晶系之外的情况,VolLetter为'A'或'I'时,如果空间群标记符号第一个或第二个字符为空，则直接将tsgn赋值给show变量
               || VolLetter == 'I')
      {
        if (   tsgn->Extension[0] == '\0'
            || tsgn->Extension[1] == '\0') // 如果第一个或第二个字符为空，则直接将tsgn赋值给show变量
          show = tsgn;
      }

/*#######################################  4444444444444444444444  #######################################################*/
      else  // 其他情况则直接将tsgn赋值给show变量
        show = tsgn;
/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 开始根据show和show_later的情况执行  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
      if (show)  // 如果show地址不为空，则
      {
        putc(' ', fpout); // 以标准输出流打印空格
        PrintTabSgNameEntry(show, 1, 0, fpout); // 以标准输出流打印空间群标记符号，完整空间群以及HallSymbol等
        putc('\n', fpout); // 以标准输出流打印换行符
        show = NULL;

        if (show_later) // 如果show_later地址不为空则继续以标准输出流打印
        {
          putc(' ', fpout); // 如果show_later地址不为空
          PrintTabSgNameEntry(show_later, 1, 0, fpout); // 打印空间群标记符号，完整空间群以及HallSymbol等
          putc('\n', fpout);
          show_later = NULL;
        }
      }
    }
  }
}


static void ListCIF(FILE *fpout)                               // 负责CIF相关的文件写入, FILE is a special pointer.
{
  int                n;
  const char         **loop, *lbl;
  const T_TabSgName  *tsgn;                                    // TabSgName

/*字符串数组, 记录了单斜晶系相关的空间群符号内容*/
  static const char *loop_monoclinic_extensions[] =            // used for what ?
    {
  "_monoclinic_extension   # cf. _symmetry_space_group_id",
  "_monoclinic_axis        # cf. IT Vol. A 1983 sec. 2.16.",
  "_monoclinic_setting     # cf. IT Vol. A 1983 tab. 2.16.1.",
  "_monoclinic_cellchoice  # cf. IT Vol. A 1983 sec. 2.16.(i) & fig. 2.6.4.",
  "",
  " b   b  abc   1",
  " b1  b  abc   1",
  " b2  b  abc   2",
  " b3  b  abc   3",
  "-b   b  c-ba  1",
  "-b1  b  c-ba  1",
  "-b2  b  c-ba  2",
  "-b3  b  c-ba  3",
  " c   c  abc   1",
  " c1  c  abc   1",
  " c2  c  abc   2",
  " c3  c  abc   3",
  "-c   c  ba-c  1",
  "-c1  c  ba-c  1",
  "-c2  c  ba-c  2",
  "-c3  c  ba-c  3",
  " a   a  abc   1",
  " a1  a  abc   1",
  " a2  a  abc   2",
  " a3  a  abc   3",
  "-a   a  -acb  1",
  "-a1  a  -acb  1",
  "-a2  a  -acb  2",
  "-a3  a  -acb  3",
      NULL
    };

/*字符串数组, 记录了所有晶系相关的一些内容*/
  static const char *loop_symmetry_space_group[] =
    {
  "_symmetry_space_group_id",
  "_symmetry_space_group_name_sch",
  "_symmetry_space_group_name_h-m   # recognised IUCr CIF data names",
  "_symmetry_space_group_name_hall  # recognised IUCr CIF data names",
      NULL
    };


  fprintf(fpout, "data_ notation\n\n");                              // print related string into file

  fprintf(fpout, "loop_\n");

  for (loop = loop_monoclinic_extensions; *loop; loop++) {           // loop the loop_monoclinic_extensions str array, and print related element into fpout.
    if ((*loop)[0]) fprintf(fpout, "    %s", *loop);                 
    putc('\n', fpout);
  }

  putc('\n', fpout);
  putc('\n', fpout);

  fprintf(fpout, "loop_\n");

  for (loop = loop_symmetry_space_group; *loop; loop++) {
    if ((*loop)[0]) fprintf(fpout, "    %s", *loop);                 // 遍历 loop_symmetry_space_group 字符串数组，以标准输出流打印所有内容
    putc('\n', fpout);
  }

  putc('\n', fpout);

  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)                   // TabSgName record all 230 space group and some related information, 遍历空间群列表
  {
    n = fprintf(fpout, "    %3d", tsgn->SgNumber);                   // print the space group number into fpout, and return the parameter numbers 

    if (tsgn->Extension[0])
      n += fprintf(fpout, ":%s", tsgn->Extension);                   // 空间群标记符号first地址不为空时写入空间群标记符号

    if (tsgn->SgNumber < 1 || tsgn->SgNumber > 230) {
      SetSgError("Internal Error: ListCIF()");                       // 当空间群数字符号超出1到230时，标记信息错误
      return;
    }

    while (n < 14) { putc(' ', fpout); n++; }                        
    putc(' ', fpout); n++;

    n += fprintf(fpout, "%s", SchoenfliesSymbols[tsgn->SgNumber]);   // 打印SchoenfliesSymbols符号 into fpout

    while (n < 22) { putc(' ', fpout); n++; }
    putc(' ', fpout); n++;

    n += PrintFullHM_SgName(tsgn, '_', fpout);                       // 打印完整空间群信息

    while (n < 36) { putc(' ', fpout); n++; }
    putc(' ', fpout);

    for (lbl = tsgn->HallSymbol; *lbl; lbl++)                        // HallSymbol
    {
      if (*lbl == ' ' && lbl != tsgn->HallSymbol)  
        putc('_', fpout);
      else
        putc(*lbl, fpout);
    }

    putc('\n', fpout);
  }
}


static void PutAllXYZ(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  const char    *xyz;
  char          buf0[8], buf1[8], buf2[8];


  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      if (nLoopInv > 1 || nTrV > 1)
        putc('#', fpout);

      if (nTrV > 1)
        fprintf(fpout, " +(%s %s %s)",
          FormatFraction(TrV[0], STBF, 0, buf0, sizeof buf0 / sizeof (*buf0)),
          FormatFraction(TrV[1], STBF, 0, buf1, sizeof buf1 / sizeof (*buf1)),
          FormatFraction(TrV[2], STBF, 0, buf2, sizeof buf2 / sizeof (*buf2)));

      if (nLoopInv > 1)
        fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

      if (nLoopInv > 1 || nTrV > 1)
        putc('\n', fpout);

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        for (i = 0; i < 9; i++)
          SMx.s.R[i] =              f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

            xyz = RTMx2XYZ(&SMx, 1, STBF, 0, 0, 1, ", ", NULL, 0);
        if (xyz)
          fprintf(fpout, "%s\n", xyz);
        else
        {
          SetSgError("Internal Error: PutAllXYZ()");
          return;
        }
      }
    }
  }

  putc('\n', fpout);
}


static void PutMaple(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  int           iMatrix;
  char          buf0[8], buf1[8], buf2[8];


  iMatrix = 0;

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      if (nLoopInv > 1 || nTrV > 1)
        putc('#', fpout);

      if (nTrV > 1)
        fprintf(fpout, " +(%s %s %s)",
          FormatFraction(TrV[0], STBF, 0, buf0, sizeof buf0 / sizeof (*buf0)),
          FormatFraction(TrV[1], STBF, 0, buf1, sizeof buf1 / sizeof (*buf1)),
          FormatFraction(TrV[2], STBF, 0, buf2, sizeof buf2 / sizeof (*buf2)));

      if (nLoopInv > 1)
        fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

      if (nLoopInv > 1 || nTrV > 1)
        putc('\n', fpout);

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        for (i = 0; i < 9; i++)
          SMx.s.R[i] =              f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

        fprintf(fpout, "m%d", ++iMatrix);
        PrintMapleRTMx(&SMx, 1, STBF, NULL, fpout);
      }
    }
  }

  putc('\n', fpout);
}


static void PutSpaceSymFile(const T_SgInfo *SgInfo, FILE *fpout)
{
  unsigned int       SgID;
  int                iList, SuppressMx, f, i;
  int                nTrV, iTrV, nLoopInv, iLoopInv;
  const int          *TrV;
  const T_RTMx       *lsmx;
  const T_TabSgName  *tsgn;


      tsgn = SgInfo->TabSgName;
  if (tsgn && tsgn->SgLabels == NULL) tsgn = NULL;

  SgID = 0;

  if (tsgn != NULL)
    SgID = SgID_Number(tsgn);

  fprintf(fpout, "%u '", SgID);

  if (tsgn != NULL)
    PrintFullHM_SgName(tsgn, 0, fpout);
  else if (SgInfo->HallSymbol[0])
    fprintf(fpout, "%s", SgInfo->HallSymbol);
  else
    fprintf(fpout, "Unknown");

  putc('\'', fpout);
  putc('\n', fpout);

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

      iList = SgInfo->OrderL;
  if (iList > 1)
  {
    iList--;
    SuppressMx = 1;
  }
  else
    SuppressMx = 0;

  fprintf(fpout, "%d\n", iList);

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        if (SuppressMx == 0)
        {
          for (i = 0; i < 3; i++)
            fprintf(fpout, " %12.8f %12.8f %12.8f %12.8f\n",
              (double) f * lsmx->s.R[3 * i + 0],
              (double) f * lsmx->s.R[3 * i + 1],
              (double) f * lsmx->s.R[3 * i + 2],
              (double) iModPositive(f * lsmx->s.T[i] + TrV[i], STBF) / STBF);

          putc(':',  fpout);
          putc('\n', fpout);
        }

        SuppressMx = 0;
      }
    }
  }
}


static void PutShelx(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           Latt_N = 0, iList;
  const T_RTMx  *lsmx;
  const char    *xyz;


  if (SgInfo->InversionOffOrigin != 0)
    fprintf(fpout, "***WARNING***: %s\n",
      "Shelx manual: the origin MUST lie on a center of symmetry");

  switch (SgInfo->LatticeInfo->Code)
  {
    case 'P': Latt_N = 1; break;
    case 'A': Latt_N = 5; break;
    case 'B': Latt_N = 6; break;
    case 'C': Latt_N = 7; break;
    case 'I': Latt_N = 2; break;
    case 'R':
              if (SgInfo->ExtraInfo == EI_Obverse)
                Latt_N = 3; break;
    case 'S':
    case 'T':
              SetSgError("Shelx supports R-obverse only");
              return;
    case 'F': Latt_N = 4; break;
    default:
      goto ReturnError;
  }

  /* N must be made negative if the structure is non-centrosymmetric
   */
  if (SgInfo->Centric != -1)
    Latt_N = -Latt_N;

  fprintf(fpout, "LATT %2d\n", Latt_N);

  lsmx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
        xyz = RTMx2XYZ(lsmx, 1, STBF, 1, 1, 0, ", ", NULL, 0);
    if (xyz)
      fprintf(fpout, "SYMM %s\n", xyz);
    else
      goto ReturnError;
  }

  putc('\n', fpout);

  return;

  ReturnError:

  SetSgError("Internal Error: PutShelx()");
  return;
}


static void PutSchakal(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           iList, nMx, i;
  int           nTrV, iTrV;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  const char    *xyz;


  if (Sg_nLoopInv(SgInfo) == 2)
    fprintf(fpout, "DU -x,-y,-z\n");

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  if (nTrV > 1)
  {
    fprintf(fpout, "DU");

    InitRotMx(SMx.s.R, 1);

    TrV += 3;

    for (iTrV = 1; iTrV < nTrV; iTrV++, TrV += 3)
    {
      for (i = 0; i < 3; i++)
        SMx.s.T[i] = TrV[i];

          xyz = RTMx2XYZ(&SMx, 1, STBF, 0, 0, 1, ",", NULL, 0);
      if (xyz)
      {
        if (iTrV > 1)
          fprintf(fpout, " ;");

        fprintf(fpout, " %s", xyz);
      }
      else
      {
        putc('\n', fpout);
        goto ReturnError;
      }
    }

    putc('\n', fpout);
  }

  nMx = 0;

  lsmx = &SgInfo->ListSeitzMx[1];

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
        xyz = RTMx2XYZ(lsmx, 1, STBF, 0, 0, 1, ",", NULL, 0);
    if (xyz)
    {
      if (nMx % 4 == 0)
      {
        if (nMx) putc('\n', fpout);
        fprintf(fpout, "SY %s", xyz);
      }
      else
        fprintf(fpout, " ; %s", xyz);
    }
    else
    {
      putc('\n', fpout);
      goto ReturnError;
    }

    nMx++;
  }

  if (nMx)
    putc('\n', fpout);

  putc('\n', fpout);

  return;

  ReturnError:

  SetSgError("Internal Error: PutSchakal()");
  return;
}


/*此函数用于打印简单的hkl列表，参数为初始化之后的SgInfo结构体指针, hkl的相应的最大值*/
static void Simple_hklList(T_SgInfo *SgInfo, int Maxh, int Maxk, int Maxl,
                           int ListSysAbsent)   // ListSysAbsent 的实参为 F_Verbose = 1
{
  int        h, k, l, iList, restriction, M, n, i;
  int        Minh, Mink, Minl;
  int        uvw[3];   // 定义int型数组
  int        CCMx_PL[9], deterCCMx_LP = 0, hP, kP, lP;  // 定义数组以四个 int 型变量，并给其中一个赋值


  if (SgInfo->LatticeInfo->Code != 'P')  // 查看是否为P格子，如果不是则需要确定格子的对称中心
  {  // deterRotMx() returns the determinant of RotMx.
    deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP); // CCMx_LP is a pointer to a constant integer array holding the "change-of-centering" (rotation) matrix from the centred to the primitive setting. 
                 InverseRotMx(SgInfo->CCMx_LP, CCMx_PL); // deterRotMx(): returns the determinant of RotMx.
// InverseRotMx() computes the inverse of the RotMx 3x3 rotation matrix, multiplied by the determinant of RotMx, and stores the result in InvRotMx.
    if (deterCCMx_LP < 1)
      goto ReturnError;
  }

  SetListMin_hkl(SgInfo, Maxk, Maxl, &Minh, &Mink, &Minl);  // set the Min hkl of hkl list

  fprintf(stdout, ">Begin hklList\n"); // 标准输出到磁盘

  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(SgInfo, h, k, l, &restriction); // IsSysAbsent_hkl() tests if the reflection with indices hkl is systematic absent for the space group as described by SgInfo. 
                                                           // If so, a code for the symmetry operation which causes the systematic absence is returned. 
    if (SgError != NULL)
      return;

    M = BuildEq_hkl(SgInfo, NULL, h, k, l); // BuildEq_hkl() generates the symmetry equivalent indices (for the space group as described by SgInfo) of reflection hkl and stores them in the Eq_hkl structure.
                                            // Of the Friedel mates only one is included in the list.
    if (SgError != NULL)
      return;

    if (iList == 0)
    { 
      // Inside a loop over Minh|k|l...Maxh|k|l IsSuppressed_hkl() decides, whether a symmetry equivalent reflex of hkl has already appeared before. 
      // If so, a code for the corresponding symmetry operation is returned; otherwise 0 is returned.
      if ((iList = IsSuppressed_hkl(SgInfo, Minh, Mink, Minl,
                                                  Maxk, Maxl,
                                               h,    k,    l)) != 0)
        n = fprintf(stdout, "# %3d %3d %3d  %3d  [%d]",   // when hkl indics already appeared, then stout output
                            h, k, l, M, iList); 
      else
        n = fprintf(stdout, "  %3d %3d %3d  %3d",       // if it appeared first time, then stout print 
                            h, k, l, M);

      if (restriction >= 0)
      {
        while (n < 27) { n++; putc(' ', stdout); }         // print blank when n < 27
        n += fprintf(stdout, " %2d/%d", restriction, STBF);  // if TH_Restriction is not the NULL pointer, and the reflection hkl has no phase restriction, *TH_Restriction is set to -1. 
                                                             // For values >= 0 the restricted phase angles in degrees are given by P = (*TH_Restriction) * (180 / STBF), and P + 180, respectively.
      }

      while (n < 34) { n++; putc(' ', stdout); }   // when n < 34, print blank
      if (Is_si(SgInfo, h, k, l) == 1)  // Is_si() returns 1 if the phase angle of the reflex with index hkl is a structure semi-invariant, 0 otherwise. Set_si() has to be called before the first call to Is_si().
        n += fprintf(stdout, " s.i.");

      while (n < 41) { n++; putc(' ', stdout); }
      Set_uvw(SgInfo, h, k, l, uvw);
      for (i = 0; i < SgInfo->n_si_Vector; i++)
        n += fprintf(stdout, " %3d", uvw[i]);

      if (SgInfo->LatticeInfo->Code != 'P')
      {
        hP = h * CCMx_PL[0] + k * CCMx_PL[3] + l * CCMx_PL[6];
        kP = h * CCMx_PL[1] + k * CCMx_PL[4] + l * CCMx_PL[7];
        lP = h * CCMx_PL[2] + k * CCMx_PL[5] + l * CCMx_PL[8];

        if (hP % deterCCMx_LP || kP % deterCCMx_LP || lP % deterCCMx_LP)
          goto ReturnError;

        hP /= deterCCMx_LP;
        kP /= deterCCMx_LP;
        lP /= deterCCMx_LP;

        while (n < 55) { n++; putc(' ', stdout); }
          n += fprintf(stdout, " P  %3d %3d %3d",
                               hP, kP, lP);
      }

      putc('\n', stdout);
    }
    else if (ListSysAbsent)
      fprintf(stdout, "# %3d %3d %3d  %3d  (%d)\n",
                      h, k, l, M, iList);
  }

  fprintf(stdout, ">End hklList\n");

  return;

  ReturnError:

  SetSgError("Internal Error: Simple_hklList()");
  return;
}


/* ****************************************************************************
   some code for the handling of lattice constants
 */


#include <math.h>


typedef struct {
                 double   a, b, c;
                 double   alpha, beta, gamma;
                 double   sa, sb, sg;
                 double   ca, cb, cg;
                 double   v;
                 char     calcs, calcc;
               }
               T_LatticeConstants;


#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif

#define PIover180 (M_PI / 180.)

#define EpsPI (1.e-6) /* ARBITRARY */


static double sinC(double arg)
{
  if (M_PI_2 - EpsPI <= arg && arg <= M_PI_2 + EpsPI)
    return 1.;
  else
    return sin(arg);
}


static double cosC(double arg)
{
  if (M_PI_2 - EpsPI <= arg && arg <= M_PI_2 + EpsPI)
    return 0.;
  else
    return cos(arg);
}


static int Lc2RLc(T_LatticeConstants *lc, T_LatticeConstants *rlc)
{
  /* Transformation Lattice Constants -> Reciprocal Lattice Constants
     after Kleber, W., 17. Aufl., Verlag Technik GmbH Berlin 1990, P.352
   */

  double  D;


  if (lc->calcs)
  { lc->sa = sinC(lc->alpha); lc->sb = sinC(lc->beta); lc->sg = sinC(lc->gamma);
    lc->calcs = 0;
  }

  if (lc->calcc)
  { lc->ca = cosC(lc->alpha); lc->cb = cosC(lc->beta); lc->cg = cosC(lc->gamma);
    lc->calcc = 0;
  }

  D = 1. - lc->ca * lc->ca - lc->cb * lc->cb - lc->cg * lc->cg
         + 2. * lc->ca * lc->cb * lc->cg;
  if (D < 0.) return -1;

  lc->v = lc->a * lc->b * lc->c * sqrt(D);
  if (lc->v == 0.) return -1;

  if (lc->sa == 0. || lc->sb == 0. || lc->sg == 0.) return -1;

  if (rlc != NULL)
  {
    rlc->a = lc->b * lc->c * lc->sa / lc->v;
    rlc->b = lc->c * lc->a * lc->sb / lc->v;
    rlc->c = lc->a * lc->b * lc->sg / lc->v;
    rlc->ca = (lc->cb * lc->cg - lc->ca) / (lc->sb * lc->sg);
    rlc->cb = (lc->cg * lc->ca - lc->cb) / (lc->sg * lc->sa);
    rlc->cg = (lc->ca * lc->cb - lc->cg) / (lc->sa * lc->sb);
    rlc->alpha = acos(rlc->ca);
    rlc->beta  = acos(rlc->cb);
    rlc->gamma = acos(rlc->cg);
    rlc->sa = sinC(rlc->alpha);
    rlc->sb = sinC(rlc->beta);
    rlc->sg = sinC(rlc->gamma);
    rlc->v = 1. / lc->v;
    rlc->calcs = 0;
    rlc->calcc = 0;
  }

  return 0;
}


static void Lc2MetricalMx(T_LatticeConstants *lc, double *G)
{
  G[0] =        lc->a * lc->a;
  G[1] = G[3] = lc->a * lc->b * lc->cg;
  G[2] = G[6] = lc->a * lc->c * lc->cb;

  G[4] =        lc->b * lc->b;
  G[5] = G[7] = lc->b * lc->c * lc->ca;

  G[8] =        lc->c * lc->c;
}


static int HarmonizeSgLatCon(T_SgInfo *SgInfo, T_LatticeConstants *lc, int np)
{
  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      if (np != 6) goto IllUnitCell;
      break;
    case XS_Monoclinic:
      if (np != 4 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': if (np != 6) lc->beta  = lc->alpha;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->gamma = lc->alpha;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Orthorhombic:
      if (np != 3 && np != 6) goto IllUnitCell;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Tetragonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b; break;
        case 'y': lc->c = lc->a; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a; break;
        default:
          goto IntErr;
      }
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Trigonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      if (SgInfo->UniqueDirCode == '*')
      {
        if (np != 6) lc->alpha = lc->b * PIover180;
        lc->c = lc->b = lc->a;
        lc->gamma = lc->beta = lc->alpha;
        break;
      }
    case XS_Hexagonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b;
                  lc->alpha = 120. * PIover180;
                  lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': lc->c = lc->a;
                  lc->beta  = 120. * PIover180;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a;
                  lc->gamma = 120. * PIover180;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Cubic:
      if (np != 1 && np != 6) goto IllUnitCell;
      lc->c = lc->b = lc->a;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    default:
      goto IntErr;
  }

  return  0;

  IntErr: SetSgError("Internal Error: HarmonizeSgLatCon()");
  return -1;

  IllUnitCell: SetSgError("Error: Illegal UnitCell or SpaceGroup");
  return -1;
}


static void MxMultiply(double *ab, double *a, double *b, int ma, int na, int nb)
{
  int     i, j, k;
  double  *ai, *aij, *bk, *bkj;

  ai = a;

  for (i = 0; i < ma; i++)
  {
    bk = b;

    for (k = 0; k < nb; k++)
    {
      aij = ai;
      bkj = bk;

      *ab = 0.;

      for (j = 0; j < na; j++)
      {
        *ab += (*aij) * (*bkj);

        aij++;
        bkj += nb;
      }

      ab++;
      bk++;
    }

    ai += na;
  }
}


static int TransformLatticeConstants(T_LatticeConstants *LatConA,
                                     int np,
                                     T_LatticeConstants *LatConB,
                                     T_SgInfo *SgInfo,
                                     int *InvCBMxR)
{
  int     i, j;
  double  GA[9], GB[9], GAR[9], R[9], Rt[9];


  if (HarmonizeSgLatCon(SgInfo, LatConA, np) != 0)
    return -1;

  LatConA->calcs = 1;
  LatConA->calcc = 1;

  /* just to check LatConA and to compute sin and cos of angles
   */
  if (Lc2RLc(LatConA, LatConB) != 0) {
    SetSgError("Error: Illegal UnitCell");
    return -1;
  }

  Lc2MetricalMx(LatConA, GA);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
       R[i * 3 + j] = InvCBMxR[i * 3 + j] / (double) CRBF;
      Rt[i * 3 + j] = InvCBMxR[j * 3 + i] / (double) CRBF;
    }

  MxMultiply(GAR, GA, R, 3, 3, 3);
  MxMultiply(GB, Rt, GAR, 3, 3, 3);

  if (GB[0] < 0. || GB[4] < 0. || GB[8] < 0.)
    goto ReturnError;

  LatConB->a = sqrt(GB[0]);
  LatConB->b = sqrt(GB[4]);
  LatConB->c = sqrt(GB[8]);

  LatConB->alpha = GB[5] / LatConB->b / LatConB->c;
  LatConB->beta  = GB[2] / LatConB->c / LatConB->a;
  LatConB->gamma = GB[1] / LatConB->a / LatConB->b;

  if (   LatConB->alpha < -1. || LatConB->alpha > 1.
      || LatConB->beta  < -1. || LatConB->beta  > 1.
      || LatConB->gamma < -1. || LatConB->gamma > 1.)
    goto ReturnError;

  LatConB->alpha = acos(LatConB->alpha);
  LatConB->beta  = acos(LatConB->beta );
  LatConB->gamma = acos(LatConB->gamma);

  LatConB->calcs = 1;
  LatConB->calcc = 1;

  return 0;

  ReturnError:

  SetSgError("InternalError: Corrupt InvCBMxR");
  return -1;
}


/* ****************************************************************************
 */


static void usage(void)  // 在命令行参数解析错误的情况下返回usage()内容
{
  static const char *quick_help[] =
    {
      "-Hall|VolA|VolI   select conventions",
      "-ListTable[=#]    print [parts of] internal table",
      "-CIF              print internal table in CIF format",
      "-XYZ              print something like \"-x, y+1/2, z\"",
      "-AllXYZ           print all symmetry operations",
      "-Maple            print symmetry matrices in Maple format",
      "-Space            print symmetry file for AVS SpaceModule",
      "-Shelx            print Shelx LATT & SYMM cards",
      "-Schakal          print Schakal DU & SY cards",
      "-hklList          print simple hkl listing",
      "-Standard         compute transformation to \"standard\" setting",
    "-UnitCell=\"a..g\"  unit cell constants a, b, c, alpha, beta, gamma",
      "-v                be more verbose",
      "-Verify           debug option: verify transformations",
      "-ClearError       debug option: clear errors and continue",
      NULL
    };

  const char  **qh;


  fprintf(stderr,
    "usage: %s [options] [SpaceGroupName_or_# [SpaceGroupName_or_#]]\n",
    progn);

  for (qh = quick_help; *qh; qh++)
    fprintf(stderr, "  %s\n", *qh);

  putc('\n', stderr);

  fprintf(stderr, "examples: %s 68\n",                          progn);
  fprintf(stderr, "          %s C2/m:c2 -XYZ\n",                progn);
  fprintf(stderr, "          %s \"Oh^3\" -Shelx\n",             progn);
  fprintf(stderr, "          %s -Hall \"-F 4y 2\" -Standard\n", progn);
  fprintf(stderr, "          %s -VolI 15 -VolA 15\n",           progn);
  fprintf(stderr, "          %s -ListTable=68\n",               progn);

  exit(1);
}


static void ShowCBMx(T_RTMx *CBMx, T_RTMx *InvCBMx, int F_Maple)
{
  if (F_Maple) {
    PrintMapleRTMx(   CBMx, CRBF, CTBF, "   CBMx", stdout);
    PrintMapleRTMx(InvCBMx, CRBF, CTBF, "InvCBMx", stdout);
  }
  else {
    fprintf(stdout, "   CBMx = %s\n",
      RTMx2XYZ(   CBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
    fprintf(stdout, "InvCBMx = %s\n",
      RTMx2XYZ(InvCBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
  }
}


typedef struct
  {
    int                Convention;
    const char         *SgName;
    const T_TabSgName  *InpTSgN;
    const T_TabSgName  *RefTSgN;
    T_RTMx             CBMx, InvCBMx;                        // Rotation-Translation matrics
  }
  T_SgList;                                                  // T_SgList struct 


int main(int argc, char *argv[])                             // main function
{
  int                 i, n, HaveSpace, pos_hsym;             // HaveSpace, used for situation that there is space in command line parameters; pos_hsym, result that return by ParseHallSymbol(SgName, SgInfo)
  int                 F_Convention, Last_F_Convention;       // the type (convention) of space group
  int                 F_ListTable, F_CIF;                    // parts of internal table, internal table in CIF format
  int                 F_XYZ, F_AllXYZ, F_Maple;              // symmetry operation, all symmetry operations, symmetry matrices in Maple format
  int                 F_Space, F_Shelx, F_Schakal;
  int                 F_hklList;                             // simple hkl listing
  int                 F_Standard, F_UnitCell;                // standard transformation, unit cell constants
  int                 F_Verbose, F_Verify, F_ClearError;     // be more verbose, verify transformations
  T_LatticeConstants  LatConA, LatConB;                      // unit cell constant
  char                *cp, xtrac; 
  const char          *SgName;                               // define the struct of spg name
  const T_TabSgName   *tsgn;                                 // T_TabSgName struct, stores the information of Hallsymbol, spg number, extenstion, sglabel
  T_SgInfo            SpgrInfo[2], BC_SgInfo, *SgInfo;       // T_SgInfo struct, stores most information of spg
  int                 nSgList, iSgList;                      // work with T_SgList
  T_SgList             SgList[2];                            // stores Convention, SgName, T_TabSgNmae, T_RTMx
  T_RTMx              *CBMx, *InvCBMx;                       // change-of-Basis-matrices
  T_RTMx              CCBMx, CInvCBMx;


/*
  Macintosh extras (Courtesy Jon Tischler <TischlerJZ@ornl.gov>)
 */
#ifdef __THINK__
  console_options.nrows = CONSOLE_LINES;
  console_options.ncols = CONSOLE_COLUMNS;
  console_options.title = "\psgInfo   version 1.0.1";
#endif
#ifdef __MWERKS__
  SIOUXSettings.autocloseonquit = FALSE;
  SIOUXSettings.asktosaveonclose = TRUE;
  SIOUXSettings.columns = CONSOLE_COLUMNS;
  SIOUXSettings.rows = CONSOLE_LINES;
#endif
#if defined(__THINK__) || defined(__MWERKS__)
  argc = ccommand(&argv);
#endif

/*init parameters*/
  nSgList = 0;                                               // nubmer of sglist

  F_Convention = 'A'; Last_F_Convention = 0;
  F_ListTable = 0;
  F_CIF = 0;
  F_XYZ = 0;
  F_AllXYZ = 0;
  F_Maple = 0;                                               // symmetry matrices in Maple format
  F_Space = 0;
  F_Shelx = 0;
  F_Schakal = 0;
  F_hklList = 0;
  F_Standard = 0;
  F_UnitCell = 0;
  F_Verbose = 0;
  F_Verify = 0;
  F_ClearError = 0;

/* loop command line paras  */
// check for the usage
  for (i = 1; i < argc; i++)                                // str_icmp() function used to check whether two str are equal
  {
    if      (str_icmp(argv[i], "-Hall") == 0) {             // check the convention of space group
      F_Convention = 'H';
      Last_F_Convention = 0;
    }
    else if (str_icmp(argv[i], "-VolA") == 0) {
      F_Convention = 'A';                                    // VolA
      Last_F_Convention = 'A';
    }                                                                                           
    else if (   str_icmp(argv[i], "-VolI") == 0
    
             || str_icmp(argv[i], "-Vol1") == 0) {  
      F_Convention = 'I';                                    // VolI
      Last_F_Convention = 'I';
    }

    /* check -ListTable, print */
    else if (str_ibegin(argv[i], "-ListTable") == 0)   // to check whether this element is matched with '-ListTable'
    {
                cp = argv[i] + 10;                     // adress，10 may repersent number of '-ListTable'
      if      (*cp == '\0') 
        F_ListTable = -1;                              // when there is nothing behand '-ListTable' print all number of space、HallSymbol、and table of other two symbol
      else if (*cp++ == '=')                           // print info of specific space group 
      {
        n = sscanf(cp, "%d%c", &F_ListTable, &xtrac);  // extract spg number 
                                                            
        if (n != 1 || F_ListTable <   1
                   || F_ListTable > 230) usage();      // validate the spg number
      }
      else               // print usage for situation expect '-ListTable' and '-ListTable=xxx'
        usage();
    }

// check whether there are related parameters
    else if (str_icmp(argv[i], "-CIF") == 0)                 // print internal table in CIF format               
      F_CIF = 1;

    else if (str_icmp(argv[i], "-XYZ") == 0)                 // print something like "-x, y+1/2, z"
      F_XYZ = 1;

    else if (str_icmp(argv[i], "-AllXYZ") == 0)              // print all symmetry operations
      F_AllXYZ = 1;

    else if (str_icmp(argv[i], "-Maple") == 0)               // print symmetry matrices in Maple format
      F_Maple = 1;

    else if (str_icmp(argv[i], "-Space") == 0)               // print symmetry file for AVS SpaceModule
      F_Space = 1;

    else if (str_icmp(argv[i], "-Shelx") == 0)               // print Shelx LATT & SYMM cards
      F_Shelx = 1;

    else if (str_icmp(argv[i], "-Schakal") == 0)             // print Schakal DU & SY cards
      F_Schakal = 1;

    else if (str_icmp(argv[i], "-hklList") == 0)             // print simple hkl listing
      F_hklList = 1;

    else if (str_icmp(argv[i], "-Standard") == 0)            // compute transformation to "standard" setting
      F_Standard = 1;

    else if (str_ibegin(argv[i], "-UnitCell=") == 0)           // sscanf the unit parameters behand '-UnitCell'
    {                                                          // unit cell constants a, b, c, alpha, beta, gamma
      F_UnitCell = sscanf(&argv[i][10], "%lf%lf%lf%lf%lf%lf", 
        &LatConA.a,     &LatConA.b,    &LatConA.c,
        &LatConA.alpha, &LatConA.beta, &LatConA.gamma);        // sscanf unit cell with lf format

      if (F_UnitCell < 1)     // when unit cell recording returns error
        usage();

      if (F_UnitCell > 3) LatConA.alpha *= PIover180; 
      if (F_UnitCell > 4) LatConA.beta  *= PIover180;
      if (F_UnitCell > 5) LatConA.gamma *= PIover180;
    }
    else if (str_icmp(argv[i], "-v") == 0)                     // verbose means to show the detail, be more verbose 
      F_Verbose = 1;

    else if (str_icmp(argv[i], "-Verify") == 0)                // debug option: verfy transformations
      F_Verify = 1;

    else if (str_icmp(argv[i], "-ClearError") == 0)            // debug option: clear errors and continue
      F_ClearError = 1;

/*the space group information must defined in the first two command line paramters*/
    else if (nSgList < 2)    // check for space group when other options do not match
    {
      SgName = argv[i];      // check for spg name

      while (*SgName == ' ' || *SgName == '\t') SgName++;      // remove ' ' and '\t' in parameter

      if (F_Convention == 'H' && isdigit(*SgName))             // check again for the convention
      /*SgList: T_SgList, stores Convention, T_TabSgName, T_RTMX*/
        SgList[nSgList].Convention = 'A';   // set the convention in the sglist
      else
        SgList[nSgList].Convention = F_Convention;            

      SgList[nSgList].SgName  = SgName;                        // set the SgName
      SgList[nSgList].InpTSgN = NULL;                          // init the input name and reference name of spg
      SgList[nSgList].RefTSgN = NULL;

      nSgList++;
    }
    else
      usage();                                                 // nSgList>=2, should be only 2 argv paras
  }

/*print all or part of internal table based on requirement*/
  if (F_ListTable) 
  { /*ListTabSgName() function print symbols of spg*/
    ListTabSgName(F_ListTable, Last_F_Convention, stdout);     // F_ListTable should be-1 or between 1 and 230
    PrintClearSgError(1, 0);                                   // check whether there is any error ? PrintClearSgError(int ClearError, int CertainSgError)
    putc('\n', stdout);
  }

  if (F_CIF) 
  {
    ListCIF(stdout);                                           // print list into CIF 
    PrintClearSgError(1, 0);                                   
    putc('\n', stdout);
  }

  if (nSgList == 0)                                            // can not locate the space group
  {
    if (F_ListTable == 0 && F_CIF == 0)                        // and can not find command line parameters of F_ListTable and F_CIF
      usage();
    else                                                       // when one of -ListTable and -CIF or both of them are defined, just exit 
      exit(0);                                                 // check me later !!!
  }

  if (F_Space == 0)                                            // print symmetry file for AVS SpaceModule when F_Space = 1
  {
    putc('#', stdout);                                         // print '#' in first line

    for (i = 0; i < argc; i++)
    {
      putc(' ', stdout);

      HaveSpace = 0;

      if (i) {
        for (n = 0; argv[i][n]; n++) {                         // check whether there is space in ith commmand line parameter
          if (isspace(argv[i][n])) {                           // isspace() check wheter str is ' ' 
            HaveSpace = 1;                                     // break when ' ' exist
            break;
          }
        }
      }

      if (HaveSpace == 0)                                      // print to stdout when there is no ' ' in ith argv parameters 
        fprintf(stdout, "%s", argv[i]);
      else 
      {                                                        // if there is space in ith command line parameter, then print all elements in stdout.
        putc('"', stdout);

        for (n = 0; argv[i][n]; n++)
          if (argv[i][n] == '"') putc('+',        stdout);
          else                   putc(argv[i][n], stdout);

        putc('"', stdout);
      }
    }

    putc('\n', stdout);
  }

  BC_SgInfo.MaxList = 0;                                       // set the MaxList to 0, usually use 192
  BC_SgInfo.ListSeitzMx = NULL;                                // init ListSeitzMx 
  BC_SgInfo.ListRotMxInfo = NULL;                              // init ListRotMxInfo

  for (iSgList = 0; iSgList < nSgList; iSgList++)              // loop SgList
  {
    if (iSgList) putc('\n', stdout);

    if (nSgList > 1 || F_Standard)                             // when there are two choices
      fprintf(stdout, "Setting %c:\n\n", "AB"[iSgList]);

    SgInfo = &SpgrInfo[iSgList];                               // SpgrInfo is struct array of T_SgInfo

    SgInfo->MaxList = 192;                                     // set the MaxList in SgInfo

    SgInfo->ListSeitzMx
      = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));   // allocate memory for ListSeitzMx 
    if (SgInfo->ListSeitzMx == NULL) NotEnoughCore();

#ifndef No_ListRotMxInfo                                           // check me later !!!
    SgInfo->ListRotMxInfo
      = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo)); // allocate memory for list rotation matrix
    if (SgInfo->ListRotMxInfo == NULL) NotEnoughCore();
#else
    SgInfo->ListRotMxInfo = NULL;
#endif

/*SgList is the struct array of T_SgList, T_SgList struct contains Convention, T_TabSgNmae and T_RTMx struct*/
// parse the parameters in SgList
    F_Convention = SgList[iSgList].Convention;                   // parse convention, such as 'A' and 'I'
    SgName          = SgList[iSgList].SgName;                    // parse spg name

    tsgn = NULL;                                                 // int tsgn struct (T_TabSgName) 
    // tsgn stores symbols of a spg including Hall and others

    if (F_Convention == 'A' || F_Convention == 'I')              // when F_convention is 'A' or 'I'
    {
          tsgn = FindTabSgNameEntry(SgName, F_Convention);       // try to find HallSymbol
      if (tsgn == NULL)
      {
        PrintClearSgError(1, 0);
        progerror("Error: Unknown Space Group Symbol");
      }

      if (F_Space == 0)                                          // print spg info into stdout
      {
        fprintf(stdout, "Space Group  ");
        PrintTabSgNameEntry(tsgn, 0, 0, stdout);                 // print space group related information
        putc('\n', stdout);
      }

      SgName = tsgn->HallSymbol;                                 // find HallSymbol, and assign it to SgName
    }

    SgList[iSgList].InpTSgN = tsgn;                              // InpTSgN is T_TabSgName struct. 

    InitSgInfo(SgInfo);                                          // init SgInfo Struct

    SgInfo->TabSgName = tsgn;                                    //  assign tsgn to SgInfo
    if (tsgn) SgInfo->GenOption = 1;                             // if tsgn is not NULL, asign SgInfo->GenOption as 1 ??? check me later !!!

    pos_hsym = ParseHallSymbol(SgName, SgInfo);                  // parse HallSymbol, to find all possible symmetry operations

    if (SgError != NULL)                                         // check SgError 
    {
      fprintf(stdout, "    %s\n", SgName);
      for (i = 0; i < pos_hsym; i++) putc('-', stdout);          // puts result and error
      fprintf(stdout, "---^\n");
      fprintf(stdout, "%s\n", SgError);
      exit(1);
    }

    if (CompleteSgInfo(SgInfo) != 0)                             // complete SgInfo()
      PrintClearSgError(F_ClearError, 1);

    if (tsgn == NULL && F_Space == 0)                            // tsgn is NULL，and F_Space is 0
    {
      if (SgInfo->TabSgName)                                     // check SgInfo->TabSgName, print tabagname info into stdout
      {
        fprintf(stdout, "Space Group  ");
        PrintTabSgNameEntry(SgInfo->TabSgName, 0, 0, stdout);
        putc('\n', stdout);
      }
      else
        fprintf(stdout, "Hall Symbol  %s\n", SgInfo->HallSymbol); // otherwise HallSymbol
    }

    PrintClearSgError(F_ClearError, 0);                           // check error, if F_ClearError if is not 0

#if USE_GS_SI
    if (Try_GS_si(SgInfo) < 0)
#else
    if (Set_si(SgInfo) < 0)
#endif
      PrintClearSgError(F_ClearError, 1);

    if (F_Space == 0) {
      ListSgInfo(SgInfo, F_XYZ, F_Verbose, stdout);               // if F_Space is0，print info spg into stdout
      PrintClearSgError(F_ClearError, 0);                         // which means information similar to international table
    }

    if (F_AllXYZ) {
      PutAllXYZ(SgInfo, stdout);                                  // print all symmetry operations
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Maple) {                                                // print symmetry matrics in Maple format into stdout
      PutMaple(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Space) {                                                // print symmetry filr for AVS SpaceModule into stdout
      PutSpaceSymFile(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Shelx) {                                                // print Shelx LATT & SYMM cards into stdout
      PutShelx(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Schakal) {                                              // print Schakal DU & SY cards into stdout
      PutSchakal(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_hklList) {                                              // print simplp hkl listing into stdout
      Simple_hklList(SgInfo, 4, 4, 4, F_Verbose);
      PrintClearSgError(F_ClearError, 0);
    }

    if (nSgList > 1 || F_Standard)                                // when nSgList>1 or F_Standard is not 0, means there is another setting
    {                                                             // or to get the standard setting
         CBMx = &SgList[iSgList].CBMx;                            // record change-of-basis matrices; check me later !!!
      InvCBMx = &SgList[iSgList].InvCBMx;

      SgList[iSgList].RefTSgN = FindReferenceSpaceGroup(SgInfo,
                                                        CBMx, InvCBMx);   // find the reference TabSgName with CBMx and InvCBMx
      PrintClearSgError(F_ClearError, 0);

      if (SgList[iSgList].RefTSgN)                                        // if reference space group was found
      {
        if (F_Verbose || F_Verify)                                        // more details are required or need to verify the transformation
        {
          fprintf(stdout, "Change of Basis => Reference Setting  ");      // print the table of sgname entry of reference space group
          PrintTabSgNameEntry(SgList[iSgList].RefTSgN, 0, 0, stdout);
          putc('\n', stdout);

          ShowCBMx(CBMx, InvCBMx, F_Maple);                               // and show the CBMx and InvCBMx
          PrintClearSgError(F_ClearError, 0);
        }

        if (F_Verify)                                                     // if verify is required
        {
          if (BC_SgInfo.MaxList == 0)
          {
            BC_SgInfo.MaxList = 192;

            BC_SgInfo.ListSeitzMx                                                  // malloc memory for BC_SgInfo
              = malloc(BC_SgInfo.MaxList * sizeof (*BC_SgInfo.ListSeitzMx));
            if (BC_SgInfo.ListSeitzMx == NULL) NotEnoughCore();

            BC_SgInfo.ListRotMxInfo
              = malloc(BC_SgInfo.MaxList * sizeof (*BC_SgInfo.ListRotMxInfo));
            if (BC_SgInfo.ListRotMxInfo == NULL) NotEnoughCore();
          }

          InitSgInfo(&BC_SgInfo);                                                   // init the sginfo of BC_SgInfo
// TODO: check CBMx and InvCBMx later !!!
          if (TransformSgInfo(SgInfo, CBMx, InvCBMx, &BC_SgInfo) == 0)              // call the TransformSgInfo() function to check whether two spg is same
            CompleteSgInfo(&BC_SgInfo);                                             // complete the BC_SgInfo

          if (SgError)                                                              // print error if there is error
          {
            PrintClearSgError(F_ClearError, 0);
            SgList[iSgList].RefTSgN = NULL;
          }
          else if (BC_SgInfo.TabSgName != SgList[iSgList].RefTSgN)                  // check the BC_sginfo whether can match with the reference sginfo 
          {
            fprintf(stdout, "Hall Symbol  %s\n", BC_SgInfo.HallSymbol);
            SetSgError("Verify Error: Wrong CBMx/InvCBMx");
            PrintClearSgError(F_ClearError, 0);
            SgList[iSgList].RefTSgN = NULL;
          }
          else
            fprintf(stdout, "Verify O.K.\n\n");                                      // if matched, just print ok
        }
      }

          tsgn = SgList[iSgList].RefTSgN;                                            // get the reference tsgn
      if (tsgn && F_Standard && nSgList == 1)                                        // change SgList[i] to standard setting
      {
        if (Last_F_Convention == 'A' || Last_F_Convention == 'I')
          SgList[nSgList].Convention = Last_F_Convention;
        else
          SgList[nSgList].Convention = 'A';

        SgList[nSgList].SgName  = SchoenfliesSymbols[tsgn->SgNumber];
        SgList[nSgList].InpTSgN = NULL;
        SgList[nSgList].RefTSgN = NULL;

        nSgList++;                                                                   // ??? check me later !!!
      }
    }
  }

  if (   nSgList == 2
      && SgList[0].RefTSgN &&           SgList[1].RefTSgN
      && SgList[0].RefTSgN->SgNumber == SgList[1].RefTSgN->SgNumber)                  // when two settings matched with each other 
  {
    putc('\n', stdout);
    fprintf(stdout, "Change of Basis Setting A -> Setting B:\n");                        

    RTMxMultiply(   &CCBMx, &SgList[1].InvCBMx, &SgList[0].CBMx,
                 CRBF, CRBF * CTBF);                                                  // ???
    RTMxMultiply(&CInvCBMx, &SgList[0].InvCBMx, &SgList[1].CBMx,
                 CRBF, CRBF * CTBF);

    for (i = 0; i < 12; i++)
    {
      if (   CCBMx.a[i] % CRBF) break;                                                // check me later !!!
      if (CInvCBMx.a[i] % CRBF) break;

         CCBMx.a[i] /= CRBF;
      CInvCBMx.a[i] /= CRBF;
    }

    if (i < 12)
    {
      SetSgError("Internal Error: Can't combine CBMx's");
      PrintClearSgError(1, 1);
    }

    else
    {
      ShowCBMx(&CCBMx, &CInvCBMx, F_Maple);
      PrintClearSgError(F_ClearError, 0);

      if (F_Verify)
      {
        InitSgInfo(&BC_SgInfo);

        if (TransformSgInfo(&SpgrInfo[0], &CCBMx, &CInvCBMx, &BC_SgInfo) == 0)
          CompleteSgInfo(&BC_SgInfo);

        if (SgError)
          PrintClearSgError(F_ClearError, 1);

        else if (strcmp(SpgrInfo[1].HallSymbol, BC_SgInfo.HallSymbol) != 0)
        {
          fprintf(stdout, "Hall Symbol  %s\n", SpgrInfo[1].HallSymbol);
          fprintf(stdout, "Hall Symbol  %s\n", BC_SgInfo.HallSymbol);
          SetSgError("Verify Error: Wrong CBMx/InvCBMx");
          PrintClearSgError(F_ClearError, 1);
        }
        else
          fprintf(stdout, "Verify O.K.\n");
      }

      if (F_UnitCell)
      {
        putc('\n', stdout);

        if (TransformLatticeConstants(&LatConA, F_UnitCell,
                                      &LatConB, &SpgrInfo[0],
                                      CInvCBMx.s.R) != 0)
          PrintClearSgError(0, 1);

        fprintf(stdout,
          "Setting A UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
          LatConA.a, LatConA.b, LatConA.c,
          LatConA.alpha / PIover180,
          LatConA.beta  / PIover180,
          LatConA.gamma / PIover180);

        fprintf(stdout,
          "Setting B UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
          LatConB.a, LatConB.b, LatConB.c,
          LatConB.alpha / PIover180,
          LatConB.beta  / PIover180,
          LatConB.gamma / PIover180);
      }
    }

    putc('\n', stdout);
  }

  exit(0); /* old VAX didn't like "return 0;" */
  return 0;
}
