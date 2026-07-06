"""Python package for the PnetCDF commux C ABI.

The package distributes the public C headers under ``pinc/include`` and
``libpnc`` under ``pinc/lib``. The optional Python API is a thin binding over
the same C ABI used by C callers.
"""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    __version__ = version("pinc")
except PackageNotFoundError:
    __version__ = "0.0.0"

package_root = Path(__file__).resolve().parent
include_dir = package_root / "include"
library_dir = package_root / "lib"

try:
    from ._pnetcdf_commux import (
        ReduceOp,
        allreduce_int64,
        barrier,
        finalize,
        init,
        init_env,
        rank,
        recv_int64,
        send_int64,
        size,
    )
except ImportError:
    pass

__all__ = [
    "include_dir",
    "library_dir",
    "ReduceOp",
    "init",
    "init_env",
    "finalize",
    "rank",
    "size",
    "barrier",
    "allreduce_int64",
    "send_int64",
    "recv_int64",
]
