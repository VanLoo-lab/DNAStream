import os
import pathlib
import time
import functools

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





def require_file_exists(func):
    """Decorator to ensure input file exists before calling the function."""
    @functools.wraps(func)
    def wrapper(self, fname, *args, **kwargs):
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File '{fname}' does not exist.")
        return func(self, fname, *args, **kwargs)
    return wrapper