import os, sys, time, platform
print("python", sys.version.split()[0], platform.machine())
print("cpus", os.cpu_count())
try:
    import psutil
    vm = psutil.virtual_memory()
    print("RAM total GB", round(vm.total/2**30,1), "avail GB", round(vm.available/2**30,1))
except Exception as e:
    print("psutil n/a", e)
    import ctypes
    class MS(ctypes.Structure):
        _fields_=[("dwLength",ctypes.c_ulong),("dwMemoryLoad",ctypes.c_ulong),("ullTotalPhys",ctypes.c_ulonglong),("ullAvailPhys",ctypes.c_ulonglong),("ullTotalPageFile",ctypes.c_ulonglong),("ullAvailPageFile",ctypes.c_ulonglong),("ullTotalVirtual",ctypes.c_ulonglong),("ullAvailVirtual",ctypes.c_ulonglong),("ullAvailExtendedVirtual",ctypes.c_ulonglong)]
    m=MS(); m.dwLength=ctypes.sizeof(MS); ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    print("RAM total GB", round(m.ullTotalPhys/2**30,1), "avail GB", round(m.ullAvailPhys/2**30,1))
print("zlib", __import__("zlib").ZLIB_VERSION)
try:
    import numpy; print("numpy", numpy.__version__)
except Exception as e: print("numpy", e)
