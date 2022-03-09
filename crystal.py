import numpy as np
import numpy.ctypeslib as npct
import ctypes as ct
from sginfo_table import PG_Names, XS_Name, EI_Name, SchoenfliesSymbols, STBF, CRBF, CTBF, SRBF
import pandas as pd
sginfo = ct.CDLL("./libsginfo.so")
import hkl_IO as hkl_io
import time
import os
import plot

# UnitCell struct
class UnitCell(ct.Structure):
    _fields_ = [("a", ct.c_double),
                ("b", ct.c_double),
                ("c", ct.c_double),
                ("alpha", ct.c_double),
                ("beta", ct.c_double),
                ("gamma", ct.c_double)]

def init_unitcell(a, b, c, alpha, beta, gamma):
    """to init unit cell parameters into UnitCell struct."""
    a, b, c = ct.c_double(a), ct.c_double(b), ct.c_double(c)
    alpha, beta, gamma = ct.c_double(alpha), ct.c_double(beta), ct.c_double(gamma)
    return UnitCell(a, b, c, alpha, beta, gamma)

def I1overI2_scale2file(fp1, fp2, fp2_scale, unitcell):
    """function to scale hkl2 against hkl1, and return a scaled file of hkl2"""
    sginfo.scale_two_file_I1overI2.argtypes = [ct.c_char_p, ct.c_char_p, ct.c_char_p, UnitCell]
    sginfo.scale_two_file_I1overI2.restype = ct.c_int
    fp1, fp2, fp2_scale = ct.c_char_p(fp1.encode('utf-8')), ct.c_char_p(fp2.encode('utf-8')), ct.c_char_p(fp2_scale.encode('utf-8'))
    return sginfo.scale_two_file_I1overI2(fp1, fp2, fp2_scale, unitcell)

def scale_merge_multi_files(filepath, unitcell, datatype='3D-EDT'):
    """to load multi data from filepath and scale them,
    then to merge them all"""
    if datatype == '3D-EDT':
        fns = hkl_io.load_hkl(filepath, ext='.hkl')
    elif datatype == 'XDS_SHELX':
        fns = hkl_io.load_hkl(filepath, ext='.HKL')
    elif datatype == 'CAP':
        pass
    elif datatype == 'PETS2':
        pass
    elif datatype == 'GENERAL':
        fns = hkl_io.load_hkl(filepath, ext='.hkl')
    print(fns)

    fns_ge = []                                        # to store hkl data file name with 'general' format
    fns_scale = []                                     # stores scaled files
    if os.path.exists(fns[0]):
        parent_dir = os.path.dirname(fns[0])
        path_general = parent_dir + '/general'
        if not os.path.exists(path_general):
            os.makedirs(path_general)
        path_scale = parent_dir + '/scale'
        if not os.path.exists(path_scale):
            os.makedirs(path_scale)

        for i in range(len(fns)):                      # loops all data and write them all into general format
            if os.path.exists(fns[i]):
                basename = os.path.basename(fns[i])
                root, ext = os.path.splitext(basename)
                data = hkl_io.read_hkl(fns[i], datatype=datatype)
                print(f"Successed to load file: {fns[i]}.")
                fn_general = path_general + '/' + root + '_general.hkl'
                fn_scal = path_scale + '/' + root + '_scaled.hkl'
                fns_ge.append(fn_general)
                fns_scale.append(fn_scal)
                hkl_io.produce_hkl(data, fn= fn_general, datatype='general')
            else:
                print(f"file {fns[i]} was failed to load, check you data!")
                return
        if len(fns_ge) == len(fns) == len(fns_scale):
            for i in range(len(fns_ge)):
                if i == 0:
                    data = hkl_io.read_hkl(fns_ge[i], datatype='GENERAL')           # write the reference data into general format.
                    hkl_io.produce_hkl(data, fn= fns_scale[0], datatype='general')
                elif os.path.exists(fns_ge[i]):
                    fp1 = fns_ge[0]
                    fp2 = fns_ge[i]
                    fp2_scale = fns_scale[i]
                    I1overI2_scale2file(fp1, fp2, fp2_scale, unitcell)
            path_out = path_scale + '/Merged.hkl'
            return hkl_io.df_mean_value_merge(fns_scale, out=path_out)
        else:
            print("Error for scaling files! Check your data again!")
    else:
        print("hkl files was failed to load!")
        return

def gene_unique_reflection_list(spg, filename_unique='unique.txt', maxh=35, maxk=35, maxl=35):
    """To generate independent reflection list.
    spg: space group name, e.g. 'Pnma' or '62'
    Return: a hkl list recording independent reflections with a limit of maxh, maxk, maxl.
    """
    sginfo.hkl_list.argtypes = (ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int)
    sginfo.hkl_list.restype = ct.c_int
    spg = ct.c_char_p(spg.encode('utf-8'))
    filename_unique = ct.c_char_p(filename_unique.encode('utf-8'))
    h, k, l = ct.c_int(maxh), ct.c_int(maxk), ct.c_int(maxl)
    return sginfo.hkl_list(spg, filename_unique, maxh, maxk, maxl)

def find_hkl_repeted(fn, fn_out, spg):
    sginfo.hkl2unique.argtypes = (ct.c_char_p, ct.c_char_p, ct.c_char_p)
    sginfo.hkl2unique.restype = ct.c_int
    spg = ct.c_char_p(spg.encode('utf-8'))
    fn = ct.c_char_p(fn.encode('utf-8'))
    fn_out = ct.c_char_p(fn_out.encode('utf-8'))
    return sginfo.hkl2unique(fn,  fn_out, spg)

class T_LatticeInfo(ct.Structure):
    """
    Code :    'P','A','B','C','I','R','S','T','F'
    nTrVector: 1,  2,  2,  2,  2,  3,  3,  3,  4
    """
    _fields_ = [("Code", ct.c_char),
                ("nTrVector", ct.c_int),
                ("TrVector", ct.POINTER(ct.c_int))]

class s(ct.Structure):
    _fields_ = [("R", ct.c_int * 9),
                ("T", ct.c_int * 3)]


class T_RTMx(ct.Union):
    _fields_ = [("s", s),
                ("a", ct.c_int * 12)]


class T_RotMxInfo(ct.Structure):
    """
    EigenVector: the rotation axis of the affiliated rotation Matrix

    Order: one of {1, 2, 3, 4, 6, -1, -2, -3, -4, -6}

    Inverse: 0 - rotation matrix has a positive turn angle
            -1 - rotation matrix has a negative turn angle

    RefAxis: one of {o, x, y, z}

    Dircode: one of {. = " ' | \\ *}
    """
    _fields_ = [("EigenVector", ct.c_int * 3),
                ("Order", ct.c_int),
                ("Inverse", ct.c_int),
                ("RefAxis", ct.c_char),
                ("DirCode", ct.c_char)]


class T_TabSgName(ct.Structure):
    _fields_ = [("HallSymbol", ct.c_char_p),
                ("SgNumber", ct.c_int),
                ("Extension", ct.c_char_p),
                ("SgLabels", ct.c_char_p)]

class T_SgInfo(ct.Structure):
    _fields_ = [("GenOption", ct.c_int),
                ("Centric", ct.c_int),
                ("InversionOffOrigin", ct.c_int),
                ("LatticeInfo", ct.POINTER(T_LatticeInfo)),
                ("StatusLatticeTr", ct.c_int),
                ("OriginShift", ct.c_int * 3),
                ("nList", ct.c_int),
                ("MaxList", ct.c_int),
                ("ListSeitzMx", ct.POINTER(T_RTMx)),
                ("ListRotMxInfo", ct.POINTER(T_RotMxInfo)),
                ("OrderL", ct.c_int),
                ("OrderP", ct.c_int),
                ("XtalSystem", ct.c_int),
                ("UniqueRefAxis", ct.c_int),
                ("UniqueDirCode", ct.c_int),
                ("ExtraInfo", ct.c_int),
                ("PointGroup", ct.c_int),
                ("nGenerator", ct.c_int),
                ("Generator_iList", ct.c_int * 4),
                ("HallSymbol", ct.c_char * 40),
                ("TabSgName", ct.POINTER(T_TabSgName)),
                ("CCMx_LP", ct.POINTER(ct.c_int)),
                ("n_si_Vector", ct.c_int),
                ("si_Vector", ct.c_int * 9),
                ("si_Modulus", ct.c_int * 3)]


class T_Eq_hkl(ct.Structure):
    """
    M: Multiplicity
    N: Number of equivalent hkl to follow
    If hkl == 000, M=N=1; if hkl != 000, M=2*N.
    """
    _fields_ = [("M", ct.c_int),
                ("N", ct.c_int),
                ("h", ct.c_int * 24),
                ("k", ct.c_int * 24),
                ("l", ct.c_int * 24),
                ("TH", ct.c_int * 24)]

def check_hkl_absent(spg, h, k, l):
    """check whether hkl indices is absent reflection."""
    sginfo.check_hkl_absent.argtypes = [ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
    sginfo.check_hkl_absent.restype = ct.c_int
    spg = ct.c_char_p(spg.encode('utf-8'))
    h, k, l = ct.c_int(h), ct.c_int(k), ct.c_int(l)
    absent = sginfo.check_hkl_absent(spg, h, k, l)
    return absent

def are_sym_equivalent_hkl(SgInfo, h1, k1, l1, h2, k2, l2):
    """check whether two hkl indices are symmetry equivalent reflections."""
    sginfo.AreSymEquivalent_hkl.argtypes = [ct.POINTER(T_SgInfo), ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int]
    sginfo.AreSymEquivalent_hkl.restype = ct.c_int
    h1, k1, l1 = ct.c_int(h1), ct.c_int(k1), ct.c_int(l1)
    h2, k2, l2 = ct.c_int(h2), ct.c_int(k2), ct.c_int(l2)
    return sginfo.AreSymEquivalent_hkl(SgInfo, h1, k1, l1, h2, k2, l2)

def get_sginfo(spg):
    """get T_SgInfo struct from a spg symbol with defalut VolLetter A setting."""
    sginfo.get_sginfo.argtypes = [ct.c_char_p]
    sginfo.get_sginfo.restype = T_SgInfo
    spg = ct.c_char_p(spg.encode('utf-8'))
    return sginfo.get_sginfo(spg)

def get_tsgn(spg):
    """get T_TabSgName struct from a spg symbol with defalut VolLetter A setting."""
    sginfo.get_tsgn.argtypes = [ct.c_char_p]
    sginfo.get_tsgn.restype = ct.POINTER(T_TabSgName)
    spg = ct.c_char_p(spg.encode('utf-8'))
    return sginfo.get_tsgn(spg)

def get_Eq_hkl(spg, h, k, l):
    """get T_Eq_hkl struct from a spg and hkl indices with defalut VolLetter A setting."""
    sginfo.get_Eq_hkl.argtypes = [ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
    sginfo.get_Eq_hkl.restype = T_Eq_hkl
    spg = ct.c_char_p(spg.encode('utf-8'))
    h, k, l = ct.c_int(h), ct.c_int(k), ct.c_int(l)
    return sginfo.get_Eq_hkl(spg, h, k, l)

def PG_Index(sginfo):
    """return PG_Index from T_SgInfo struct."""
    pg_code = sginfo.PointGroup
    return pg_code // (33*12)

def point_group(sginfo):
    """return point group name from T_SgInfo struct."""
    pg_index = PG_Index(sginfo)
    return PG_Names[pg_index]

def get_laue(spg):
    """return laue group from space group."""
    sginfo.get_laue.argtypes = [ct.c_char_p]
    sginfo.get_laue.restype = ct.c_char_p
    spg = ct.c_char_p(spg.encode('utf-8'))
    return sginfo.get_laue(spg).decode('utf-8')

def unique_axis(sginfo):
    """return unique axis."""
    if sginfo.UniqueRefAxis!= 0 or sginfo.UniqueRefAxis!= ord('o'):
        return chr(sginfo.UniqueRefAxis)
    else:
        return None

def arr_merge2unique(arr):
    """to merge reflection intensity with same indices.
    must to find symmetry equivalent reflections first."""
    df = hkl_io.hkl_arr2df(arr)
    val = df.groupby(df.index).mean()
    val = val.reset_index()
    val = val.values
    val1 = np.array(list(val[:, 0])).astype(float)
    val2 = np.array(val[:, 1])
    val3 = np.array(val[:, 2])
    val = np.c_[val1, val2.T]
    val = np.c_[val, val3.T]
    return val

def arr_merge2unique_file(arr, out='Merged_unique.hkl', outtype='general_no_nomalize'):
    """to merge reflection intensity with same indices into a file.
    must to find symmetry equivalent reflections first."""
    val = arr_merge2unique(arr)
    if outtype == 'general_no_nomalize':
        outtype = 'general'
        hkl_io.produce_hkl(val, out, datatype=outtype)
    elif outtype == 'general_nomalize':
        outtype = 'general'
        arr_3max = val[:, 3].max()
        val[:, 3] = 99999.99 * val[:, 3] / arr_3max
        val[:, 4] = 99999.99 * val[:, 4] / arr_3max
        hkl_io.produce_hkl(val, out, datatype=outtype)
    if outtype == 'shelx_no_nomalize':
        outtype = 'shexl'
        hkl_io.produce_hkl(val, out, datatype=outtype)
    elif outtype == 'shelx_nomalize':
        outtype = 'shelx'
        arr_3max = val[:, 3].max()
        val[:, 3] = 99999.99 * val[:, 3] / arr_3max
        val[:, 4] = 99999.99 * val[:, 4] / arr_3max
        hkl_io.produce_hkl(val, out, datatype=outtype)

class Crystal(object):
    def __init__(self, symbol, cell=[]):
        super(Crystal, self).__init__()
        self.symbol = symbol
        self.init_sginfo(symbol)
        self.init_tsgn(symbol)
        self.cell = cell
        self.parse_cell()

    def init_sginfo(self, symbol):
        self.sginfo = sginfo = get_sginfo(symbol)
        self.point_group = point_group(sginfo)
        self.laue_group = get_laue(symbol)
        self.crsytal_system = XS_Name[sginfo.XtalSystem]
        self.unique_axis = unique_axis(sginfo)
        self.extrainfo = EI_Name[sginfo.ExtraInfo]
        self.inversion_off_origin = sginfo.InversionOffOrigin
        self.order = sginfo.OrderL
        self.order_p = sginfo.OrderP
        self.centric = sginfo.Centric

    def init_tsgn(self, symbol):
        self.tsgn = tsgn = get_tsgn(symbol)
        self.hall_symbol = tsgn.contents.HallSymbol.decode('utf-8')
        self.spg_number = tsgn.contents.SgNumber
        self.setting = self.extension = tsgn.contents.Extension.decode('utf-8')
        self.schoenflies_symbol = SchoenfliesSymbols[self.spg_number]
        self.spg_label = self.hermann_mauguin = tsgn.contents.SgLabels.decode('utf-8')
        self.qualif = None                   # TODO: revise me later.
        self.isChiral = False                # TODO: revise me later.
        self.isEnantiomorphic = False        # TODO: revise me later.

    def parse_cell(self):
        if len(self.cell) == 6:
            self.a = self.cell[0]
            self.b = self.cell[1]
            self.c = self.cell[2]
            self.alpha = self.cell[3]
            self.beta = self.cell[4]
            self.gamma = self.cell[5]
            self.cell_parsed = True
        else:
            print("Error of recording unit cell parameters!")
            self.cell_parsed = False

    def init_UnitCell_struct(self):
        if self.cell_parsed:
            a, b, c, = self.a, self.b, self.c
            alpha, beta, gamma = np.radians(self.alpha), np.radians(self.beta), np.radians(self.gamma)
            self.cell_struct = init_unitcell(a, b, c, alpha, beta, gamma)
            return self.cell_struct
        else:
            print("Error of init UnitCell struct!")

    def _Eq_hkl(self, h, k, l):
        """To get equivalent reflections of original one."""
        symbol = self.symbol
        Eq_hkl = get_Eq_hkl(symbol, h, k, l)
        multiplicity = Eq_hkl.M
        num_equivalent = Eq_hkl.N
        h = Eq_hkl.h       # list
        k = Eq_hkl.k
        l = Eq_hkl.l
        TH = Eq_hkl.TH
        sym_equi_hkl = []
        TH_list = []
        for i in range(num_equivalent):
            sym_equi_hkl.append([h[i], k[i], l[i]])
        for j in range(num_equivalent):
            TH_list.append(TH[i])
        return multiplicity, num_equivalent, sym_equi_hkl, TH_list
        # print(multiplicity, num_equivalent, sym_equi_hkl, TH_list)

    def is_hkl_absent(self, index=[1,0,0]):
        """check whether a index is absent or not"""
        symbol = self.symbol
        h, k, l = index[0], index[1], index[2]
        absent = check_hkl_absent(symbol, h, k, l)
        if absent == 0:               # for non-absent reflections
            return False
        else:
            return True               # for absent reflections

    def _getMultiplicity(self, index=[1,0,0]):
        """get the multiplicity of a specific index."""
        h, k, l = index[0], index[1], index[2]
        mul, _, _, _ = self._Eq_hkl(h, k, l)
        return mul

    def num_equi_indices(self, h, k, l):
        """get the number of equivalent refletions of index hkl."""
        _, num, _, _ = self._Eq_hkl(h, k, l)
        return num

    def sym_equi_hkl(self, h, k, l, Fridel=True):
        """get the equivalent hkl index list."""
        _, _, equi, _ = self._Eq_hkl(h, k, l)
        if Fridel:
            _, _, equi_fridel, _ = self._Eq_hkl(-h, -k, -l)
            equi.extend(equi_fridel)
        return equi

    def _getPhaseRestriction(self,index=[1,0,0]):
        """get the pahse restriction list of equivalent reflecions of index hkl"""
        h, k, l = index[0], index[1], index[2]
        _, _, _, th = self._Eq_hkl(h, k, l)
        return th

    def are_sym_equivalent_hkl(self, h1, k1, l1, h2, k2, l2):
        """check whether two index are equivalent reflections or not."""
        SgInfo = self.sginfo
        result = are_sym_equivalent_hkl(SgInfo, h1, k1, l1, h2, k2, l2)
        if result:
            return True
        else:
            return False

    def _apply_along_index(self, arr, func):               # check me later !!!
        """Expects tuple/list of 3 elements or iterable

        return 1/0"""
        if isinstance(arr, pd.Index):
            return arr.map(func)
        elif isinstance(arr, pd.DataFrame):
            return arr.index.map(func)
        elif len(arr[0]) == 3:
            # assume nested list
            return np.array(map(func, arr))
        else:
            return func(arr)

    def d_value(self, index=[1, 0, 0]):
        if self.cell_parsed:
            a, b, c, = self.a, self.b, self.c
            alpha, beta, gamma = self.alpha, self.beta, self.gamma
            h, k, l = index[0], index[1], index[2]
            d = hkl_io.dhkl_cal(a, b, c, alpha, beta, gamma, h, k , l)
            return d
        else:
            print("Parse the unit cell first!")

    def arr_d_value(self, index):
        return self._apply_along_index(index, self.d_value)

    def filter_with_d_value(self, index, dmin=1.0, dmax=50.0):
        d = self.arr_d_value(index)
        d_filter = (d > dmin) & (d < dmax)
        return index[d_filter]

    def is_absent(self, index):
        """method to check whether a index is systematically absent or not."""
        return self._apply_along_index(index, self.is_hkl_absent)

    def multiplicity(self, index):
        """method to check the multiplicity of an index."""
        return self._apply_along_index(index, self._getMultiplicity)

    def phaserestriction(self, index):
        """method to check the phase restriction of an index."""
        return self._apply_along_index(index, self._getPhaseRestriction)

    def isCentrosymmetric(self):
        return self.centric == -1

    def _symmetry_operations(self, centre_check=True, nSMx_check=True):
        """Method to generate all symmetry operations."""
        sginfo = self.sginfo
        if centre_check:
            centric = sginfo.Centric
        else:
            centric = 1
        if centric == -1:
            inversion = [1, -1]
        else:
            inversion = [1]
        n_latt_vector = sginfo.LatticeInfo.contents.nTrVector
        latt_vector = sginfo.LatticeInfo.contents.TrVector
        if nSMx_check:
            nSMx = sginfo.nList
        else:
            nSMx = 1
        print("nSMX", nSMx)
        for i in range(n_latt_vector):
            for j in range(len(inversion)):
                for k in range(nSMx):
                    s = sginfo.ListSeitzMx[k].s
                    latt_center = [latt_vector[3*i+0], latt_vector[3*i+1], latt_vector[3*i+2]]
                    MX = get_seitz_MX(s, inversion[j], latt_center)
                    yield SymOp(MX)


    def gene_uniuqe_hkl_list(self, unit_cell=None, d_min=1.0):
        if unit_cell:
            unit_cell = unit_cell
        elif self.cell:
            unit_cell = self.cell
        else:
            print("Please input unit cell parameters first.")
            return
        a, b, c, alpha, beta, gamma = unit_cell[0], \
            unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5]
        spg = self.symbol
        filename_unique = 'unique.txt'
        maxh, maxk, maxl = int(a/d_min )+1, int(b/d_min )+1, int(c/d_min )+1
        gene_unique_reflection_list(spg, filename_unique=filename_unique,
                                    maxh=maxh, maxk=maxk, maxl=maxl)
        unique = np.loadtxt(filename_unique)
        return unique

    def merge_hkl2unique_arr(self, fn, datatype='3D-EDT', dmin=0.8):
        """method to merge reflections into unique numpy array with C code."""
        data = hkl_io.read_hkl(fn, datatype=datatype)
        df = hkl_io.hkl_arr2df(data)
        d_filter = self.filter_with_d_value(df, dmin=dmin)
        absent_filter = self.filter_systematic_absences(d_filter)
        absent_filter = hkl_io.hkl_df2arr(absent_filter)
        filepath = "to_find_repeated.hkl"
        hkl_io.produce_hkl(absent_filter, fn= filepath, datatype='general')
        repeated = self.find_hkl_repeted(filepath)
        unique = arr_merge2unique(repeated)
        return unique

    def merge_multi_hkl2unique_hkl(self, fp, datatype='3D-EDT', dmin=0.8):
        """merge multi reflections into unique numpy array with C code."""
        if datatype == '3D-EDT':
            fns = hkl_io.load_hkl(fp, ext='.hkl')
        elif datatype == 'XDS_SHELX':
            fns = hkl_io.load_hkl(fp, ext='.HKL')
        elif datatype == 'CAP':
            pass
        elif datatype == 'PETS2':
            pass
        print("Starting to record data.")

        fns_unique = []
        if os.path.exists(fns[0]):
            parent_dir = os.path.dirname(fns[0])
            path_unique = parent_dir + '/unique'
            if not os.path.exists(path_unique):
                os.makedirs(path_unique)

            for i in range(len(fns)):
                if os.path.exists(fns[i]):
                    basename = os.path.basename(fns[i])
                    root, ext = os.path.splitext(basename)
                    unique = self.merge_hkl2unique_arr(fns[i], datatype=datatype, dmin=dmin)
                    fn_unique = path_unique + '/' + root + '_unique.hkl'
                    fns_unique.append(fn_unique)
                    hkl_io.produce_hkl(unique, fn= fn_unique, datatype='general')

                    print("Merging data finished:", fns[i])
                else:
                    print("Error for producing unique files! Check your data !")
            return fns_unique

        else:
            print("hkl files was failed to load!")
            return

    def filter_systematic_absences(self, df):
        """Takes a reflection list and filters reflections that are absent"""
        try:
            index = df.index
        except AttributeError:
            index = df
        sel = self.is_absent(df)
        sel = (sel == False)
        return df[sel]

    def find_hkl_repeted(self, fn):
        """method to find the repeated symmetry-equivalent reflections in a hkl file,
        Return: numpy array."""
        if os.path.exists(fn):
            spg = self.symbol
            fn_out = 'test_repeated.hkl'
            find_hkl_repeted(fn, fn_out, spg)
            data = np.loadtxt(fn_out)
            data = np.array(data, dtype='float64')
            return data
        else:
            print("Error of finding repeated hkl with c routine!")

    def scale_merge_multi_files(self, filepath, datatype='3D-EDT'):
        """method to scale multi hkl data and merge them into one hkl file."""
        if self.cell_parsed:
            unitcell = self.init_UnitCell_struct()
            return  scale_merge_multi_files(filepath, unitcell, datatype=datatype)
        else:
            print("Init unit cell parameters first!")

    def completeness(self, fn, datatype='3D-EDT', dmin=1.0):
        data = hkl_io.read_hkl(fn, datatype=datatype)
        df = hkl_io.hkl_arr2df(data)
        d_filter = self.filter_with_d_value(df, dmin=dmin)
        d_filter = hkl_io.hkl_df2arr(d_filter)
        unique_cal = self.gene_uniuqe_hkl_list(d_min=dmin)
        unique_cal_df = hkl_io.arr_unique_cal2df(unique_cal)
        unique_cal_d_filter = self.filter_with_d_value(unique_cal_df, dmin=dmin)
        print((f"Theretical reflecion numbers: {len(unique_cal_d_filter)} with resolution of {dmin} A"))
        filepath = "find_repeated.hkl"
        hkl_io.produce_hkl(d_filter, fn= filepath, datatype='general')
        repeated = self.find_hkl_repeted(filepath)
        unique = arr_merge2unique(repeated)
        print(f"Unique reflecion numbers: {len(unique)} with resolution of {dmin} A")
        completeness = len(unique) / len(unique_cal_d_filter)
        print("completeness of data is:", completeness)
        return completeness

    def R_int(self, fn, datatype='3D-EDT', dmin=1.0):
        data = hkl_io.read_hkl(fn, datatype=datatype)
        df = hkl_io.hkl_arr2df(data)
        d_filter = self.filter_with_d_value(df, dmin=dmin)
        d_filter = hkl_io.hkl_df2arr(d_filter)
        filepath = "find_repeated.hkl"
        hkl_io.produce_hkl(d_filter, fn= filepath, datatype='general')
        repeated = self.find_hkl_repeted(filepath)
        repeated_df = hkl_io.hkl_arr2df(repeated)

        def r_int(x):
            res=  sum(abs(x- (x).mean())) / sum(x)
            return res

        repeated_df_new = repeated_df.copy()
        R_int = repeated_df_new.groupby(repeated_df_new.index)
        unique_rint = R_int.agg(r_int)
        rint = unique_rint['val'].mean()
        print("R_int value is:",rint)
        return rint

    def hkl_cutoff2file(self,fn, unique=False, datatype='3D-EDT', dmin=0.8, fn_cutoff='hkl_cutoff.hkl'):
        if unique:
            d_filter = self.merge_hkl2unique_arr(fn, datatype=datatype, dmin=dmin)
        else:
            data = hkl_io.read_hkl(fn, datatype=datatype)
            df = hkl_io.hkl_arr2df(data)
            d_filter = self.filter_with_d_value(df, dmin=dmin)
            d_filter = hkl_io.hkl_df2arr(d_filter)
        hkl_io.produce_hkl(d_filter, fn= fn_cutoff, datatype='shelx')

    def plot_1D_single(self,fn, unique=False, datatype='3D-EDT', dmin=0.8, wavelength=0.02508):
        if unique:
            d_filter = self.merge_hkl2unique_arr(fn, datatype=datatype, dmin=dmin)
        else:
            data = hkl_io.read_hkl(fn, datatype=datatype)
            df = hkl_io.hkl_arr2df(data)
            d_filter = self.filter_with_d_value(df, dmin=dmin)
            d_filter = hkl_io.hkl_df2arr(d_filter)
            d_value = hkl_io.hkl_arr2df(d_filter)
            d_value = self.arr_d_value(d_value)
            d_value = np.array(d_value)
            two_theta = 2*np.arcsin(wavelength / (2*d_value)) * 180 / np.pi
            print(d_filter.shape, two_theta.shape)
            plot.plot_theta_intensity(two_theta, d_filter, show_hkl=False)

    def plot_1D_compare(self,fn1, fn2, unique=False, datatype='3D-EDT', dmin=0.8, wavelength=0.02508):
        if unique:
            d_filter1 = self.merge_hkl2unique_arr(fn1, datatype=datatype, dmin=dmin)
            d_filter2 = self.merge_hkl2unique_arr(fn2, datatype=datatype, dmin=dmin)
        else:
            data1 = hkl_io.read_hkl(fn1, datatype=datatype)
            df1 = hkl_io.hkl_arr2df(data1)
            d_filter1 = self.filter_with_d_value(df1, dmin=dmin)
            d_filter1 = hkl_io.hkl_df2arr(d_filter1)
            d_value1 = hkl_io.hkl_arr2df(d_filter1)
            d_value1 = self.arr_d_value(d_value1)
            d_value1 = np.array(d_value1)
            two_theta1 = 2*np.arcsin(wavelength / (2*d_value1)) * 180 / np.pi

            data2 = hkl_io.read_hkl(fn2, datatype=datatype)
            df2 = hkl_io.hkl_arr2df(data2)
            d_filter2 = self.filter_with_d_value(df2, dmin=dmin)
            d_filter2 = hkl_io.hkl_df2arr(d_filter2)
            d_value2 = hkl_io.hkl_arr2df(d_filter2)
            d_value2 = self.arr_d_value(d_value2)
            d_value2 = np.array(d_value2)
            two_theta2 = 2*np.arcsin(wavelength / (2*d_value2)) * 180 / np.pi

            plot.plot_theta_intensity2(two_theta1, d_filter1, two_theta2, d_filter2,show_hkl=False)




if __name__ == '__main__':

    cell = [19.55, 15.01, 6.66, 90, 103.7, 90]
    sp_symbol = "P21/n"


    sp = Crystal(sp_symbol, cell)

    t0 = time.time()

    # fp = './data/mil53/mil53_lt_104/unique/*.hkl'
    # fn = './data/mil53/mil53_lt_104/data6.hkl'
    fn1 = './data/mil53/mil53_lt_104/data6.hkl'
    fn2 = './data/mil53/mil53_lt_104/data8.hkl'

    # sp.scale_merge_multi_files(fp, datatype='GENERAL')
    # sp.merge_multi_hkl2unique_hkl(fp, datatype='3D-EDT', dmin=0.8)
    # sp.R_int(fn, datatype='3D-EDT', dmin=1.0)
    # sp.completeness(fn, datatype='3D-EDT', dmin=1.0)
    # sp.plot_1D_single(fn, unique=False, datatype='3D-EDT', dmin=1.0)
    sp.plot_1D_compare(fn1, fn2, unique=False, datatype='3D-EDT', dmin=1.0, wavelength=0.02508)

    t1 = time.time()
    t = t1 - t0
    print("Total calculation time:", t)
