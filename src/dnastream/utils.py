import os
import pathlib
import time
import functools
import numpy as np
import uuid


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

    # Helper: normalize label values (handle bytes coming from HDF5)


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
