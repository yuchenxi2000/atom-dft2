import ctypes
import numpy as np
import sys
import pathlib

if sys.platform.startswith('win'):  # Windows
    libatomdft = ctypes.CDLL('libatomdft', winmode=0)
else:  # macOS, Linux
    libatomdft = np.ctypeslib.load_library(libname='libatomdft', loader_path=pathlib.Path(__file__).parent)

# int solve_radial_from_zero(double * chi, double * dchi, double E, int l, double Z, double * Vext, double r0, double rc, int N);
libatomdft.solve_radial_from_zero.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # chi
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # dchi
    ctypes.POINTER(ctypes.c_int),  # num_nodes
    ctypes.c_double,  # E
    ctypes.c_int,  # l
    ctypes.c_double,  # Z
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # Vext
    ctypes.c_double,  # r0
    ctypes.c_double,  # rc
    ctypes.c_int,  # N
]
libatomdft.solve_radial_from_zero.restypes = ctypes.c_int

# int solve_radial_from_zero_and_inf(double * chi, double * dchi, double * dE, int * num_nodes, double E, int l, double Z, double * Vext, double r0, double rc, int N);
libatomdft.solve_radial_from_zero_and_inf.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # chi
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # dchi
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # dE
    ctypes.POINTER(ctypes.c_int),  # num_nodes
    ctypes.c_double,  # E
    ctypes.c_int,  # l
    ctypes.c_double,  # Z
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # Vext
    ctypes.c_double,  # r0
    ctypes.c_double,  # rc
    ctypes.c_int,  # N
]
libatomdft.solve_radial_from_zero_and_inf.restypes = ctypes.c_int

# void solve_poisson(double * VH, double * dVH, double * density_r, double r0, double rc, int N);
libatomdft.solve_poisson.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # VH
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # dVH
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # density_r
    ctypes.c_double,  # r0
    ctypes.c_double,  # rc
    ctypes.c_int,  # N
]
libatomdft.solve_poisson.restypes = None
