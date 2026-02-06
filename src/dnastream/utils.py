from __future__ import annotations
import os
import pathlib
import time
import functools
import numpy as np
import uuid
from typing import Callable, Any
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path


def wrap_list(val):
    if type(val) is list:
        return val
    return [val]


def full_path(fname):
    return str(pathlib.Path(fname).resolve())


def timeit(func):
    """Decorator to measure execution time of a function."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()  # Start timer
        result = func(*args, **kwargs)  # Run the function
        end_time = time.perf_counter()  # End timer
        elapsed_time = end_time - start_time
        print(f"⏱ Function '{func.__name__}' took {elapsed_time:.4f} seconds")
        return result

    return wrapper


def as_str(x, *, encoding="utf-8", errors="strict") -> str:
    if isinstance(x, np.ndarray) and x.shape == ():
        x = x.item()
    if isinstance(x, np.bytes_):
        x = bytes(x)
    if isinstance(x, bytes):
        return x.decode(encoding, errors)
    if isinstance(x, np.str_):
        return str(x)
    if isinstance(x, str):
        return x
    if isinstance(x, uuid.UUID):
        return str(x)
    raise TypeError(f"Expected string-like; got {type(x).__name__}")


def norm_key(x) -> str:
    # canonical form for registry keys
    return as_str(x).strip()


def as_str_vec(arr, *, encoding="utf-8", errors="strict") -> np.ndarray:
    """
    Vectorized-ish conversion of a 1D array of string-likes to a numpy array of Python str.

    Handles:
      - fixed-length bytes dtype: kind == 'S'
      - unicode dtype: kind == 'U'
      - object dtype containing bytes/np.bytes_/str/np.str_
    """
    a = np.asarray(arr)

    # scalar -> 1-element vector
    if a.shape == ():
        a = a.reshape(1)

    if a.dtype.kind == "U":
        # already unicode
        return a.astype(str)

    if a.dtype.kind == "S":
        # fixed-length bytes -> decode in bulk
        # np.char.decode returns dtype '<U...'
        return np.char.decode(a, encoding=encoding, errors=errors)

    if a.dtype.kind == "O":
        # object array: may mix bytes and str
        # First, decode bytes-like objects only.
        out = a.astype(object, copy=True)

        # find bytes-like entries
        is_bytes = np.fromiter(
            (
                (
                    isinstance(x, (bytes, np.bytes_))
                    or (
                        isinstance(x, np.ndarray)
                        and x.shape == ()
                        and isinstance(x.item(), (bytes, np.bytes_))
                    )
                )
                for x in out
            ),
            dtype=bool,
            count=out.size,
        )

        if is_bytes.any():
            b = out[is_bytes]
            # unwrap 0-d arrays
            b = np.array(
                [
                    x.item() if isinstance(x, np.ndarray) and x.shape == () else x
                    for x in b
                ],
                dtype=object,
            )
            b = np.array(
                [bytes(x) if isinstance(x, np.bytes_) else x for x in b], dtype=object
            )
            out[is_bytes] = np.char.decode(
                np.asarray(b, dtype="S"), encoding=encoding, errors=errors
            ).astype(object)

        # now coerce everything to str (but avoid turning None into "None" if you care)
        return out.astype(str)

    # fallback: just stringify
    return a.astype(str)


def require_file_exists(func):
    """Decorator to ensure input file exists before calling the function."""

    @functools.wraps(func)
    def wrapper(self, fname, *args, **kwargs):
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File '{fname}' does not exist.")
        return func(self, fname, *args, **kwargs)

    return wrapper


def require_file_exists_static(func):
    """Decorator to ensure input file exists before calling the function."""

    @functools.wraps(func)
    def wrapper(fname, *args, **kwargs):
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File '{fname}' does not exist.")
        return func(fname, *args, **kwargs)

    return wrapper


def _qualname(fn: Callable[..., Any]) -> str:
    mod = getattr(fn, "__module__", "") or ""
    qn = getattr(fn, "__qualname__", None) or getattr(fn, "__name__", "") or ""
    return f"{mod}.{qn}" if mod else qn


def decode(x):
    return x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else x


def decode_arr(x: Any, *, encoding: str = "utf-8") -> Any:
    """Decode a single HDF5/NumPy row into pure Python types.

    - For a NumPy structured scalar (np.void) or 0-d structured array, return dict.
    - For a 1-d structured array, return list[dict].
    - Otherwise return x unchanged.

    This avoids the 'S' dtype trap where decoded strings get cast back to bytes.
    """

    def _decode_scalar_row(row: np.void) -> dict[str, object]:
        out: dict[str, object] = {}
        for name in row.dtype.names or ():
            v = row[name]
            if isinstance(v, (bytes, np.bytes_)):
                v = v.decode(encoding)
            elif isinstance(v, np.generic):
                v = v.item()
            out[name] = v
        return out

    # structured scalar
    if isinstance(x, np.void) and getattr(x.dtype, "names", None) is not None:
        return _decode_scalar_row(x)

    # structured array
    if isinstance(x, np.ndarray) and getattr(x.dtype, "names", None) is not None:
        if x.shape == ():
            return _decode_scalar_row(x[()])
        return [_decode_scalar_row(r) for r in x]

    return x


def package_version(pkg: str = "dnastream") -> str:
    try:
        return version(pkg)
    except PackageNotFoundError:
        return "0+unknown"


def resolve_path(path: str) -> str:
    return str(Path(path).resolve())


def get_file_id(path: str) -> str:
    p = Path(path).resolve()
    try:
        st = os.stat(p)
        file_id = f"{st.st_size}:{int(st.st_mtime)}"
    except Exception:
        file_id = ""
    return file_id
