import sys
import os
import csv
import getpass
import socket
import pathlib
import time
from datetime import datetime
import functools
import pandas as pd
import h5py
import numpy as np
from enum import Enum


# import json
from .index_manager import LocalIndex, GlobalIndex, DependentIndexView
from .enums import GlobalIndexName, LocalIndexName, Modalities, TreeType, SchemaGroups

from .schema import (
    SCHEMA,
    STRUCT_ARRAYS,
    LOCAL_INDEX,
    GLOBAL_INDEX,
    COPY_NUMBER_LAYER_DICT,
    get_schema_value,
)
from .datatypes import EDGE_LIST_DTYPE, ALLELE_SPECIFIC_CN_DTYPE


def wrap_list(val):
    if type(val) is list:
        return val
    return [val]


def full_path(fname):
    return str(pathlib.Path(fname).resolve())


# TODO
# - add pyclone addition
# - add patient metadata
# - LCM coordinates
# - pseudo-bulk optional layers

"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               # Short name chr:pos:ref:alt
 │   ├── data                 # Structured array: quality scores, number of callers, etc.
 │   ├── cluster              # Integer cluster assignments
 │   ├── log                  # Index modification log
 ├── sample/                  # Sample index
 │   ├── labels               # Sample names
 │   ├── data                 # Structured array: patient ID, source, location, file paths
 │   ├── cluster              # Integer cluster assignments
 │   ├── log                  # Index modification log
 ├── tree/
 │   ├── SNV_trees/  
 │   │   ├── trees            # Variable-length edge lists of clusters
 │   │   ├── data             # Likelihood, rank, method used to generate, etc.
 │   │   ├── labels           # Tree labels
 │   ├── CNA_trees/     
 │   │   ├── trees            # Variable-length edge lists (or Newick strings)
 │   │   ├── data             # Likelihood, rank, method used to generate, etc.
 │   │   ├── labels           # Tree labels
 ├── copy_numbers/
 │   ├── bulk/
 │   │   ├── labels           # Segment labels (chrom, start, end)
 │   │   ├── index            # Bulk-specific segment index
 │   │   ├── profile          # 3D array: (segment, sample, allele-specific CN)
 │   │   ├── logr             # 2D array: (segment, sample) logR values
 │   │   ├── baf              # 2D array: (segment, sample) B-allele frequency
 │   │   ├── log              # Modification log
 │   ├── lcm/
 │   │   ├── labels           # Segment labels
 │   │   ├── index            # LCM-specific segment index
 │   │   ├── profile          # 3D array: (segment, sample, allele-specific CN)
 │   │   ├── logr             # 2D array: (segment, sample) logR values
 │   │   ├── baf              # 2D array: (segment, sample) B-allele frequency
 │   │   ├── log              # Modification log
 │   ├── scdna/
 │   │   ├── labels           # Segment labels
 │   │   ├── index            # Single-cell segment index
 │   │   ├── profile          # (sample, segment) → allele CN tuple
 │   │   ├── log              # Modification log
 ├── read_counts/
 │   ├── variant              # 2D array: (SNV, sample) variant read counts
 │   ├── total                # 2D array: (SNV, sample) total read counts
 │   ├── log                  # Read count modifications log
 ├── metadata/
 │   ├── log                  # Metadata modifications log
 │   ├── sample_info          # Sample metadata
 │   ├── processing_parameters # Processing parameters used in analysis
"""


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


class DNAStream:
    """
    DNAStream is an HDF5-based data structure for efficient storage, indexing,
    and retrieval of processed DNA sequencing data across multiple sequencing modalities.

    It provides structured storage for SNV read counts, copy number profiles, and metadata
    from bulk, LCM, and single-cell sequencing. The design ensures consistency across
    different data modalities and enables efficient querying and updating.



    Instance Attributes
    -------------------
    filename : str
        Path to the HDF5 file used for data storage.
    verbose : bool
        If True, enables verbose output during operations (default: False).
    file : h5py.File
        The HDF5 file object that manages storage and retrieval.


    Methods
    -------
    __init__(filename, verbose=False, snv_dtype, sample_dtype)
        Initializes the HDF5 storage and creates necessary datasets if they do not exist.
    add_edge_list(fname, source, location=None)
        Adds read count data from a file and updates the corresponding datasets.
    add_maf_file(fname, missing_values, required_cols)
        Adds SNV information from a Mutation Annotation Format (MAF) file.
    add_maf_files(fnames, **kwargs)
        Wrapper to add multiple MAF files.
    batch_add_snvs(labels, source_file)
        Efficiently adds multiple SNVs to the index and updates the index log file.
    batch_add_samples(labels, source_file)
        Efficiently adds multiple samples to the index and updates the index log file.
    get_snv_log()
        Returns the SNV index log as a Pandas DataFrame
    get_sample_log()
        Returns the sample index log a Pandas DataFrame
    get_snv_data(indices=None)
        Retrieves structured SNV data as a Pandas DataFrame.
    get_sample_data(indices=None)
        Retrieves structured sample data as a Pandas DataFrame.
    _resize_all(m=None, n=None)
        Resizes all datasets when adding new SNVs or samples.
    close()
        Closes the HDF5 file.
    """

    def __init__(self, filename, verbose=False, id=None, sex=None, safe=True):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose
        self.safe = safe

        self.global_idx = {}
        self.local_idx = {}
        self.in_context = False

    def __enter__(self):
        """Context manager for opening the HDF5 file."""
        self.in_context = True
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.in_context = False
        self.close()

    def _connect(self, mode="a"):
        if not self.in_context:
            print(
                "⚠️ Warning: You are connecting outside of a context manager. This is not recommended as runtime errors may corrupt the file."
            )
        if not self._is_connected():
            self.file = h5py.File(
                self.filename, mode
            )  # Append mode (does not overwrite)
            if self.verbose:
                print(f"#Stream to connection {self.filename} open...")
        else:
            if self.verbose:
                print(f"#Stream to connection {self.filename} already open...")

    def connect(self):
        self._connect()

        # create index objects and bound class methods
        self._init_indices()

    def _init_indices(self):

        # Initialize built-in global indices
        for index_name in GLOBAL_INDEX:
            self.create_global_index(
                index_name.value,
                metadata_dtype=GLOBAL_INDEX[index_name]["metadata"]["dtype"],
                tracked_tables=GLOBAL_INDEX[index_name]["tracked_tables"],
            )

        for index_name in LOCAL_INDEX:
            self.create_local_index(
                index_name,
                metadata_dtype=LOCAL_INDEX[index_name]["metadata"]["dtype"],
                tracked_tables=LOCAL_INDEX[index_name]["tracked_tables"],
            )

        def is_lcm_sample(meta):
            return meta["modality"].lower() == b"lcm"

        self.copy_number_lcm_view = DependentIndexView(
            global_index=self.global_idx["sample"],
            tracked_tables=[
                ("copy_numbers/lcm/logr", 1),
                ("copy_numbers/lcm/baf", 1),
                ("copy_numbers/lcm/profile", 1),
            ],
            predicate_fn=is_lcm_sample,
            file=self.file,
            verbose=self.verbose,
        )

        def is_bulk_sample(meta):
            return meta["modality"].lower() == b"bulk"

        self.copy_number_bulk_view = DependentIndexView(
            global_index=self.global_idx["sample"],
            tracked_tables=[
                ("copy_numbers/bulk/logr", 1),
                ("copy_numbers/bulk/baf", 1),
                ("copy_numbers/bulk/profile", 1),
            ],
            predicate_fn=is_bulk_sample,
            file=self.file,
            verbose=self.verbose,
        )

        def is_scdna_sample(meta):
            return meta["modality"].lower() == b"scdna"

        self.copy_number_scdna_view = DependentIndexView(
            global_index=self.global_idx["sample"],
            tracked_tables=[
                ("copy_numbers/scdna/logr", 1),
                ("copy_numbers/scdna/baf", 1),
                ("copy_numbers/scdna/profile", 1),
            ],
            predicate_fn=is_scdna_sample,
            file=self.file,
            verbose=self.verbose,
        )

    def initialize(self, patient, sex, overwrite=False):
        """Build and connect to an HDF5 file"""
        self._connect(mode="w" if overwrite else "a")
        if patient:
            self.set_patient_id(patient)

        else:
            self.set_patient_id("")

        if sex:
            self.set_patient_sex(sex)
        else:
            self.set_patient_sex("")
        self._recursive_build(SCHEMA)
        # self._init_indices()
        self.close()

    def __str__(self):
        """To string method"""
        m = self.get_snv_size()
        n = self.get_sample_size()

        mystr = f"DNAStream object with {m} SNVs and {n} samples"
        mystr += f"\nHDF5 File: {self.filename}"
        mystr += f"\nPatient: {self.file.attrs['id']}, sex: {self.file.attrs['sex']}"

        return mystr

    def __getitem__(self, key):
        if key not in self.file:
            raise KeyError(f"Table '{key}' not found in HDF5 file.")
        return self.file[key]

    def list_tables(self):
        """List all tables in the HDF5 file."""
        return list(self.file.keys())

    def create_global_index(self, index_name, metadata_dtype=None, tracked_tables=None):
        """
        Create a new global index and bind its accessor methods.

        Parameters
        ----------
        index_name : str
            Name of the index (e.g., 'SNV', 'sample', etc.).
        metadata_dtype : numpy.dtype, optional
            The metadata data type for the index.
        tracked_tables : list of tuples, optional
            Tables and axis to be tracked with this index (table_name, axis).

        Notes
        -----
        - Adds the new index to `self.global_idx`.
        - Dynamically binds accessor methods (`get_index`, `get_clusters`, etc.).
        """

        # Initialize the GlobalIndex and store it
        self.global_idx[index_name] = GlobalIndex(
            self.file,
            f"global_index/{index_name}",
            metadata_dtype=metadata_dtype,
            tracked_tables=tracked_tables,
        )

        # Define dynamic methods to bind
        methods = {
            "get_clusters": "get_{}_clusters",
            "size": "get_{}_size",
            "get_index": "get_{}_index",
            "get_labels": "get_{}_labels",
            "add": "batch_{}_add",
            "update_clusters": "update_{}_clusters",
            "get_log": "get_{}_log",
            "get_metadata": "get_{}_metadata",
            "insert_metadata": "insert_{}_metadata",
            "label_to_idx": "{}_label_to_idx",
            "resize": "{}_resize",
            "indices_by_label": "indices_by_{}_label",
        }

        for attr, method_template in methods.items():
            method_name = method_template.format(index_name.lower())
            actual_method = getattr(self.global_idx[index_name], attr)

            @functools.wraps(actual_method)  # Preserve function signature & docstring
            def method(self, *args, actual_method=actual_method, **kwargs):
                return actual_method(*args, **kwargs)

            # Assign docstring explicitly, ensuring actual_method.__doc__ is not None
            method.__doc__ = (
                actual_method.__doc__ or ""
            ).strip() + f" (Auto-generated for {index_name})"

            # Bind the method to the instance
            setattr(self, method_name, method.__get__(self))
            # setattr(DNAStream, method_name, method)

        if self.verbose:
            print(f"#Loading global index '{index_name}' and bound methods...")

    def create_local_index(self, index_name, metadata_dtype=None, tracked_tables=None):
        """
        Create a new local index and bind its accessor methods.

        Parameters
        ----------
        index_name : str
            Name of the index (e.g., 'tree/SNV_trees', 'copy_numbers/bulk').
        metadata_dtype : numpy.dtype, optional
            The metadata data type for the index.
        tracked_tables : list of tuples, optional
            Tables and axis to be tracked with this index (table_name, axis).

        Raises
        ------
        ValueError
            If the tracked tables are not within the same group scope.

        Notes
        -----
        - Adds the new index to `self.local_idx`.
        - Dynamically binds accessor methods (`get_labels`, `get_index`, etc.).
        - Ensures that all tracked tables are in the same group scope.
        """

        # Ensure tracked tables all belong to the same base group
        if tracked_tables:
            base_groups = {table.split("/")[0] for table, _ in tracked_tables}
            if len(base_groups) > 1:
                raise ValueError(
                    f"Inconsistent tracking groups detected: {base_groups}. "
                    "All tracked tables must belong to the same base group."
                )

        # Initialize the LocalIndex and store it
        self.local_idx[index_name] = LocalIndex(
            self.file,
            f"local_index/{index_name}",
            metadata_dtype=metadata_dtype,
            tracked_tables=tracked_tables,
        )

        # Define dynamic methods to bind
        methods = {
            "get_labels": "get_{}_labels",
            "size": "get_{}_size",
            "get_index": "get_{}_index",
            "insert_metadata": "insert_{}_metadata",
            "get_metadata": "get_{}_metadata",
            "label_to_idx": "{}_label_to_idx",
            "resize": "{}_resize",
            "indices_by_label": "indices_by_{}_label",
            "add": "batch_{}_add",
        }

        for attr, method_template in methods.items():
            method_name = method_template.format(index_name.lower().replace("/", "_"))
            actual_method = getattr(self.local_idx[index_name], attr)

            @functools.wraps(actual_method)  # Preserve function signature & docstring
            def method(self, *args, actual_method=actual_method, **kwargs):
                return actual_method(*args, **kwargs)

            # Assign docstring explicitly
            method.__doc__ = (
                actual_method.__doc__ or ""
            ).strip() + f" (Auto-generated for {index_name})"

            # Bind the method to the instance
            setattr(self, method_name, method.__get__(self))

        if self.verbose:
            print(f"#Loading local index '{index_name}' and bound methods...")

    def get_dtype(self, table):
        return self[table].dtype

    def set_patient_id(self, id):
        self.file.attrs["id"] = id

    def get_patient_id(self):
        return self.file.attrs.get("id", None)

    def set_patient_sex(self, sex):
        self.file.attrs["sex"] = sex

    def get_patient_sex(self):
        return self.file.attrs.get("sex", None)

    def _recursive_build(self, schema, path=""):
        """
        Recursively traverse the schema dictionary and create datasets in the HDF5 file.

        Parameters
        ----------
        schema : dict
            The nested schema dictionary defining the dataset structure.
        path : str, optional
            The current path in the HDF5 hierarchy, used for recursive traversal.
        """
        if isinstance(schema, dict):
            for key, value in schema.items():
                if isinstance(key, Enum):
                    key_str = key.value
                else:
                    key_str = key

                new_path = (
                    f"{path}/{key_str}" if path else key_str
                )  # Construct the HDF5 path

                # if key == "index":
                #     return
                if isinstance(value, dict) and "dtype" in value:
                    # Base case: value holds dataset initialization specs
                    self.add_dataset_to_file(
                        new_path,
                        dtype=value["dtype"],
                        shape=value.get("shape", (0,)),
                        maxshape=value.get("maxshape", (None,)),
                        chunks=value.get("chunks", None),  # Apply chunking if defined
                        compression="gzip",
                        columns=(
                            list(value["dtype"].names) if key in STRUCT_ARRAYS else []
                        ),
                    )
                else:
                    # Recursive case: value is a nested dictionary, keep traversing until dtypes is among keys
                    self._recursive_build(value, new_path)

    def add_dataset_to_file(
        self,
        path,
        dtype=h5py.string_dtype("utf-8"),
        columns=[],
        source_file="",
        **kwargs,
    ):
        """
        Internal function as wrapper to create a new dataset for the file

        Parameters
        ----------

        path : str
            path to where dataset should be added
        dtype : numpy.dtype
            the simple or complex datatype of the dataset

        """
        if path not in self.file:
            if self.verbose:
                print(f"#Creating dataset {path}...")
            self.file.create_dataset(path, dtype=dtype, **kwargs)
            if columns:
                self[path].attrs["columns"] = columns
            self._log_dataset_modification(
                path, operation="create", source_file=source_file
            )
        else:
            print(f"#Warning! '{path}' dataset exists and will not be overwritten.")

    def _log_dataset_modification(self, dataset_name, operation, source_file=""):
        """
        Logs modifications to any dataset within the HDF5 file.

        Parameters
        ----------
        dataset_name : str
            Path to the dataset that was modified.
        operation : str
            The operation performed (e.g., "add", "update", "delete", "create").
        source_file : str, optional
            The source file that triggered the modification.

        Notes
        -----
        - Logs include timestamps, affected dataset, user info, and modification details.
        - Uses a dynamically growing dataset in HDF5.
        """
        log = self["metadata/log"]
        current_size = log.shape[0]
        new_size = current_size + 1
        log.resize((new_size,))

        timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S").encode("utf-8")
        user = getpass.getuser().encode("utf-8")
        hostname = socket.gethostname().encode("utf-8")
        source_file = str(pathlib.Path(source_file).resolve()).encode("utf-8")

        # Store modification log
        log[current_size] = (
            timestamp_str,
            dataset_name.encode("utf-8"),
            operation.encode("utf-8"),
            user,
            hostname,
            source_file,
        )

    def safe_mode_enable(self):
        """
        Switches DNAStream into safe mode.
        """
        self.safe = True

    def safe_mode_disable(self):
        self.safe = False

    # @timeit
    def add_read_counts(self, fname, columns=None):
        """
        Add variant and total read count data for SNVs from a file.

        This method reads a file containing variant and total read counts, add SNVs and samples to index,
        maps SNV and sample labels to indices, and updates the HDF5 dataset efficiently. Function ignores
        user-specified column names but the first four columns should contain the SNV label, sample label,
        variant read counts and total read counts.

        Parameters
        ----------
        fname : str
            Path to the input file containing read counts.
        source : str
            Sequencing modality source (e.g., "bulk", "lcm", "scdna").
        location : str, optional
            Additional sample location metadata for the samples.

        Raises
        ------
        Exception
            If an error occurs during file processing.
        """
        try:
            rc = pd.read_csv(fname)

            if columns:
                rc.rename(columns=columns, inplace=True)

            for col in ["snv", "sample", "alt", "total"]:
                if col not in rc.columns:
                    raise ValueError(f"Column '{col}' not found in the input file.")

            samples = rc["sample"].unique().tolist()

            snvs = rc["snv"].unique().tolist()
            self.batch_snv_add(snvs, source_file=fname)

            self.batch_sample_add(samples, source_file=fname)

            sample_indices_arr = np.array(
                self.indices_by_sample_label(rc["sample"]), dtype=int
            )
            snv_indices_arr = np.array(self.indices_by_snv_label(rc["snv"]), dtype=int)

            # if source:
            #     self._update_value(
            #         sample_indices,
            #         f"{GlobalIndexName.SAMPLE.value}/metadata",
            #         "source",
            #         source,
            #     )

            # if location:
            #     self._update_value(
            #         sample_indices,
            #         f"{GlobalIndexName.SAMPLE.value}/metadata",
            #         "location",
            #         location,
            #     )

            var_counts = rc["alt"].astype(int).to_numpy()
            total_counts = rc["total"].astype(int).to_numpy()

            # Update dataset in batch instead of looping
            dat = self["read_counts"]
            # unique_snv_indices = np.unique(snv_indices_arr)  # Get unique SNV indices

            unique_snv_indices = np.unique(snv_indices_arr)  # Get unique SNV indices

            for snv in unique_snv_indices:
                mask = snv_indices_arr == snv
                s_indices = sample_indices_arr[mask]
                # Ensure that the sample indices (and corresponding values) are sorted
                order = np.argsort(s_indices)
                s_indices = s_indices[order]
                var_vals = var_counts[mask][order]
                tot_vals = total_counts[mask][order]

                dat["variant"][snv, s_indices] = var_vals
                dat["total"][snv, s_indices] = tot_vals

                # if len(s_indices) == 0:
                #     continue

                # # If the sample indices are contiguous, update via slicing.
                # if np.all(np.diff(s_indices) == 1):
                #     start_idx = s_indices[0]
                #     end_idx = s_indices[-1] + 1
                #     dat["variant"][snv, start_idx:end_idx] = var_vals
                #     dat["total"][snv, start_idx:end_idx] = tot_vals
                # else:
                #     # Otherwise, update element-by-element.
                #     for si, v, t in zip(s_indices, var_vals, tot_vals):
                #         dat["variant"][snv, si] = v
                #         dat["total"][snv, si] = t

            for arr in ["variant", "total"]:
                self._log_dataset_modification(
                    f"read_counts/{arr}", operation="update", source_file=fname
                )

        except Exception:
            self.close()
            raise

    def _add_read_count(self, snv_idx, sample_idx, var=0, total=0):
        """
        Internal function to add a single read count entry. Currently not used.

        Parameters
        ----------
        source : str
            Sequencing modality (e.g., "bulk", "lcm", "scdna").
        snv_idx : int
            Index of the SNV.
        sample_idx : int
            Index of the sample.
        var : int, optional
            Variant read count (default is 0).
        total : int, optional
            Total read count (default is 0).
        """

        dat = self[f"read_counts"]
        dat["variant"][snv_idx, sample_idx] = var
        dat["total"][snv_idx, sample_idx] = total

    def get_dataset_log(self):
        """
        Retrieve the dataset log

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for SNV updates.
        """
        return self._get_data(dataset_name=f"metadata/log")

    def _get_data(self, dataset_name, indices=None):
        """
        Internal method to retrieve structured data from the HDF5 file.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to retrieve (e.g., "SNV/metadata").
        indices : list of int, optional
            The indices of the rows to retrieve. If None, retrieves all rows.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the requested dataset.

        Notes
        -----
        - The method decodes byte strings into UTF-8.
        - If `indices` is provided, only the specified rows are retrieved.
        """
        dataset = self[dataset_name]

        columns = dataset.attrs["columns"]

        if indices is not None:
            dataset = dataset[indices, :]
        else:
            dataset = dataset[:]

        df = pd.DataFrame(dataset, columns=columns)

        # Convert byte strings to UTF-8
        for col in df.select_dtypes(include=["object"]):
            if df[col].dtype == object:
                df[col] = df[col].apply(
                    lambda x: x.decode("utf-8") if isinstance(x, bytes) else x
                )

        return df

    def add_snv_data(self, indices, df, source_file=""):
        """
        Add or update SNV data in the HDF5 file.

        Parameters
        ----------
        indices : list of int
            The indices corresponding to the SNVs being added or updated.
        df : pandas.DataFrame
            A DataFrame containing the SNV data.

        Raises
        ------
        IndexError
            If indices exceed the dataset size.

        Notes
        -----
        - The dataset is resized if necessary before adding data.
        - Uses structured NumPy arrays for efficient storage.
        """
        self._add_data(
            indices,
            df,
            dataset_name=f"{GlobalIndexName.SNV.value}/metadata",
            source_file=source_file,
        )

    def add_sample_data(self, indices, df, source_file=""):
        """
        Add or update sample data in the HDF5 file.

        Parameters
        ----------
        indices : list of int
            The indices corresponding to the samples being added or updated.
        df : pandas.DataFrame
            A DataFrame containing the sample data.

        Raises
        ------
        IndexError
            If indices exceed the dataset size.

        Notes
        -----
        - The dataset is resized if necessary before adding data.
        - Uses structured NumPy arrays for efficient storage.
        """
        self._add_data(
            indices,
            df,
            dataset_name=f"{GlobalIndexName.SAMPLE.value}/metadata",
            source_file=source_file,
        )

    def _add_data(self, indices, df, dataset_name, source_file=""):
        """
        Internal method to update structured HDF5 datasets.

        Parameters
        ----------
        indices : list of int
            The indices to update in the dataset.
        df : pandas.DataFrame
            A DataFrame containing the data to insert.
        dataset_name : str
            The name of the dataset in HDF5.

        Raises
        ------
        IndexError
            If the dataset size is too small to accommodate the indices.

        Notes
        -----
        - Converts DataFrame to a structured NumPy array before updating.
        - Only modifies the specified indices to optimize performance.
        """
        dataset = self[dataset_name]
        max_index = max(indices) + 1
        if dataset.shape[0] < max_index:
            self.file.close()
            raise IndexError(
                "Invalid indices passed to add_data method, check indices and try again"
            )

        # Convert df to structured array
        structured_data = np.zeros(len(indices), dtype=dataset.dtype)
        for col in df.columns:
            structured_data[col] = df[col].to_numpy()

        # Update data at provided indices
        dataset[indices] = structured_data

        # log modifications within the internal add function
        self._log_dataset_modification(
            dataset_name, operation="update", source_file=source_file
        )

    def _update_value(self, indices, dataset_name, col, value, source_file=""):
        """
        Internal method to update a single column in an HDF5 dataset.

        Parameters
        ----------
        indices : list of int
            The indices of the rows to update.
        dataset_name : str
            The name of the dataset in HDF5.
        col : str
            The name of the column to update.
        value : array-like or scalar
            The new values to assign.

        Raises
        ------
        IndexError
            If indices exceed dataset size.

        Notes
        -----
        - If `value` is an array, it must be sorted in the same order as `indices`.
        - The affected rows are read into memory before modification.
        """
        dataset = self[dataset_name]
        max_index = max(indices) + 1
        if dataset.shape[0] < max_index:
            self.file.close()
            raise IndexError(
                "Invalid indices passed to update_data method, check indices and try again"
            )

        # Read affected rows from HDF5 into memory
        temp_data = dataset[indices]  # Read only required rows

        # Update only the specified column
        temp_data[col] = value

        # Write back only modified rows
        dataset[indices] = temp_data
        self._log_dataset_modification(
            dataset_name, operation="update", source_file=source_file
        )

    def add_maf_files(self, fnames, **kwargs):
        """
        Add multiple MAF (Mutation Annotation Format) files to the SNV index and updates the index log.

        This method iterates over a list of MAF file paths and processes each file
        using `add_maf_file`, updating the SNV index and data in the HDF5 dataset.

        Parameters
        ----------
        fnames : list of str
            A list of file paths to MAF files.
        **kwargs : dict, optional
            Additional keyword arguments to pass to `add_maf_file`.

        Notes
        -----
        - This function loops through the list of files and calls `add_maf_file` for each.
        - Handles missing values and ensures column consistency across files.
        """
        for f in fnames:
            self.add_maf_file(f, **kwargs)

    def add_maf_file(
        self,
        fname,
        missing_values=[
            "Unknown",
            "Na",
            "N/A",
            "na",
            "nan",
            "NaN",
            "NAN",
            "NONE",
            "None",
            "",
            "__UNKNOWN__",
        ],
        required_cols=[
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Entrez_Gene_Id",
        ],
    ):
        """
        Add a single MAF (Mutation Annotation Format) file to the SNV index and updates the index log.

        This function reads a MAF file, extracts SNV-related information,
        and adds it to the HDF5 dataset. If missing values or required columns
        are not present, they are handled accordingly.

        Parameters
        ----------
        fname : str
            Path to the MAF file to be processed.
        missing_values : list of str, optional
            A list of strings that should be treated as missing values (default includes common NA representations).
        required_cols : list of str, optional
            A list of required column names that should be present in the MAF file
            (default includes standard MAF mutation columns).

        Raises
        ------
        Exception
            If an error occurs during processing, the HDF5 file is closed, and an exception is raised.

        Notes
        -----
        - The function ensures that all required columns are present.
        - Any missing columns are added with `pd.NA` values.
        - SNV labels are created by concatenating the first four columns (chr:pos:ref:alt).
        - The extracted data is stored in the HDF5 file in a structured format.
        """

        # read MAF file and extract key info
        maf = pd.read_table(fname, low_memory=False)

        missing_cols = set(required_cols) - set(maf.columns)

        maf.replace(missing_values, pd.NA, inplace=True)
        if missing_cols:
            maf = maf.reindex(
                columns=maf.columns.tolist() + list(missing_cols), fill_value=pd.NA
            )

        maf["label"] = maf.iloc[:, :4].astype(str).agg(":".join, axis=1)

        column_dict = {
            "label": "label",
            "Chromosome": "chrom",
            "Start_Position": "pos",
            "End_Position": "end_pos",
            "Reference_Allele": "ref_allele",
            "Tumor_Seq_Allele2": "alt_allele",
            "Hugo_Symbol": "hugo",
            "Entrez_Gene_Id": "gene",
        }

        maf.rename(columns=column_dict, inplace=True)

        snv_labels = maf["label"].tolist()

        try:

            snv_idx = self.batch_snv_add(snv_labels, source_file=fname)

            # sort dataframe according to the newly assigned indices in DNAStream
            maf.loc[:, "snv_idx"] = maf["label"].map(snv_idx)
            maf = maf.sort_values("snv_idx")
            indices = maf["snv_idx"].tolist()
            snv_data = maf[[val for _, val in column_dict.items()]]

            self.add_snv_data(indices, snv_data, source_file=fname)

        except Exception:

            self.close()
            raise

    @staticmethod
    def _parse_file(fname, sep_word="tree", nskip=0, sep="\t"):
        """
        Parses a text file containing multiple edge lists of trees.
        @param fname: str filename to be parsed
        @param sep_word: str a word contained in the line separating distinct trees (default 'tree')
        @param nksip: int number of rows to skip before parsing
        """

        tree_list = []
        with open(fname, "r+") as file:
            new_tree = None
            for idx, line in enumerate(file):
                if idx < nskip:
                    continue
                if sep_word in line:
                    if new_tree:
                        tree_list.append(new_tree)
                    new_tree = []

                else:
                    edge = [int(e) for e in line.strip().split(sep)]
                    new_tree.append((edge[0], edge[1]))
            if new_tree:
                tree_list.append(new_tree)

        return tree_list, None

    def _expand(self, table_list, n):
        """
        Expand the size of the tables in the HDF5 file.

        Parameters
        ----------
        table_list : list of str
            List of table names to expand.
        n : int
            the amount of additional space to allocate
        """
        table_list = wrap_list(table_list)
        for table in table_list:
            old_size = self[table].shape[0]
            new_size = old_size + n
            self[table].resize((new_size,))

    def add_trees_from_edge_lists(self, tree_list, scores=None, method=""):
        """
        Manually add  list of trees to the HDF5 dataset.

        Parameters
        ----------
        tree_list : list of list of int tuples
            List of (parent, child) edges representing the tree structure.

        scores : list of float, optional
            A list of numerical scores associated with the trees (e.g., likelihood scores).
        method : str, optional
            The method used to generate the trees (e.g., "conipher"), by default "".
        source_file : str, optional
            Path to the source file from which the trees are being added (default is an empty string).
        """

        try:
            self._add_trees(tree_list, "SNV", scores, method)
        except Exception:
            self.close()
            raise

    def _add_trees(self, tree_list, tree_type, scores, method="", source_file=""):
        """
        Add a list of trees to the HDF5 file.

        Parameters
        ----------
        tree_list : list of list of int tuples
            List of (parent, child) edges representing the tree structure.
        tree_type : str
        scores : list of float, optional
            A list of numerical scores associated with the trees (e.g., likelihood scores).
        method : str, optional
            The method used to generate the trees (e.g., "conipher"), by default "".
        source_file : str, optional
            Path to the source file from which the trees are being added (default is an empty string).
        """
        table_name = f"{SchemaGroups.TREES.value}/{tree_type}"
        tree_list = wrap_list(tree_list)

        tree_dict = self.trees_snv_resize(len(tree_list), prefix="snv_tree")

        tree_indices = list(tree_dict.values())

        numtrees = len(tree_dict)
        # self._expand(table_name, numtrees)

        # Add the tree metadata to the file
        data_dtype = LOCAL_INDEX[f"trees_{tree_type}"]["metadata"]["dtype"]
        tree_dtype = self.get_dtype(f"{SchemaGroups.TREES.value}/{tree_type}")

        if not scores:
            dat = np.array(
                [(lab, method, np.nan, -1, source_file) for lab in tree_dict],
                dtype=data_dtype,
            )
        else:
            dat = np.array(
                [
                    (lab, method, scores[i], -1, source_file)
                    for i, lab in enumerate(tree_dict)
                ],
                dtype=data_dtype,
            )

        # self[f"{table_name}/metadata"][tree_indices] = dat
        # SCHEMA[DNAStream.TREE][f"{tree_type}_"]["trees"]["dtype"]
        tree_structured = np.array(
            [np.array(tree, dtype=EDGE_LIST_DTYPE) for tree in tree_list],
            dtype=tree_dtype,
        )
        self[f"{SchemaGroups.TREES.value}/{tree_type}"][tree_indices] = tree_structured

        self._log_dataset_modification(
            f"{SchemaGroups.TREES.value}/{tree_type}", "update", source_file=source_file
        )

    def _search_data_by_filename(self, dataset_name, fname):
        """
        Search for rows in a dataset where the "file" column matches a given filename.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to search within.
        fname : str
            The filename to search for.

        Returns
        -------
        list of int
            A list of indices where the filename matches the "file" column in the dataset.
            Returns an empty list if the column "file" is not present in the dataset.
        """
        df = self._get_data(dataset_name)
        if "file" in df.columns:
            return df.index[df["file"] == fname].to_list()
        return []

    def _check_safe(self, dataset_name, source_file):
        indices = self._search_data_by_filename(dataset_name, source_file)
        if len(indices) > 0:

            if self.verbose:
                print(
                    f"#Warning! Entries already exist from {source_file} in '{dataset_name}'"
                )

            return False
        else:
            return True

    def add_trees_from_file(self, fname, tree_type="SNV", method=""):
        """
        Add phylogenetic trees from a file to the HDF5 dataset.

        Parameters
        ----------
        fname : str
            Path to the file containing tree structures.
        tree_type : str, optional
            The type of tree being added (e.g., "SNV", "CNA", "clonal"), by default "SNV".
        method : str, optional
            The method used to generate the tree (e.g., "conipher"), by default "".
        safe : boolean, optional
            Safe mode to refuse add trees to data if the entries from a given filename exists, by default True


        Raises
        ------
        Exception
            If there is an error processing the file or adding trees.

        Notes
        -----
        - The function checks if trees from the source file already exist.
        - Trees are parsed based on the method (e.g., "conipher" vs. space-separated format).
        - Logs dataset modifications after adding trees.
        """
        try:
            source_file = str(pathlib.Path(fname).resolve())

            dataset_name = (
                f"local_index/{SchemaGroups.TREES.value}/{tree_type}/metadata"
            )

            # if self.safe and not self._check_safe(dataset_name, source_file):
            #     print(
            #         "#Warning! Attempting overwrite in Safe mode, use safe=F, to force append trees."
            #     )
            #     return

            # Parse tree file according to method
            if self.verbose:
                print(f"#Adding {tree_type} trees from {fname}...")
            if method == "conipher":
                tree_list, scores = self._parse_file(fname)
            else:
                tree_list, scores = self._parse_file(fname, nskip=1, sep=" ")

            self._add_trees(
                tree_list, tree_type, scores, method, source_file=source_file
            )

            if self.verbose:
                print(f"#{len(tree_list)} {tree_type} trees added from {method}.")

        except Exception:
            self.close()
            raise

    def add_snv_tree_from_edge_list(
        self, edge_list, method="", score=np.nan, rank=np.nan
    ):
        """
        Add an SNV tree from an edge list to the dataset.

        Parameters
        ----------
        edge_list : list of tuples
            List of (parent, child) edges representing the tree structure.
        method : str, optional
            The method used to generate the tree (e.g., "conipher"), by default "".
        score : float, optional
            A numerical score associated with the tree (e.g., likelihood score), by default NaN.
        rank : int, optional
            Rank of the tree among other inferred trees, by default NaN.

        Notes
        -----
        - The function logs modifications to the SNV tree dataset.
        - The tree is added with metadata including method, score, and rank.
        """
        data_dtype = LOCAL_INDEX["trees_SNV"]["metadata"]["dtype"]

        dat = np.array([("", method, score, rank, "")], dtype=data_dtype)

        self._add_trees(edge_list, "SNV", data=dat)
        self._log_dataset_modification(
            f"{SchemaGroups.TREES.value}/SNV_tree", operation="update"
        )

    def _extract_indices_by_column(self, dataset_name, name, values):
        vals = self[dataset_name][name][:]  # Load the column data
        indices = np.where(np.isin(vals, values))[0]  # Get matching indices
        return indices

    def _extract_data(self, dataset_name, snv_indices=None, sample_indices=None):
        """
        Extracts a subset of the dataset based on SNV and sample indices.

        Parameters
        ----------
        dataset_name : str
            The HDF5 dataset name.
        snv_indices : list or array-like, optional
            Indices of SNVs to extract. If None, selects all SNVs.
        sample_indices : list or array-like, optional
            Indices of samples to extract. If None, selects all samples.

        Returns
        -------
        np.ndarray
            A NumPy array containing the extracted data.
        """

        if snv_indices and sample_indices:
            arr = self[dataset_name][snv_indices, :]
            return arr[:, sample_indices]
        elif snv_indices:
            return self[dataset_name][snv_indices]
        elif sample_indices:
            return self[dataset_name][:, sample_indices]
        else:
            return self[dataset_name][:]

    def extract_read_counts(
        self,
        tables=["variant", "total"],
        sources=None,
        snv_labels=None,
        snv_indices=None,
        sample_labels=None,
        sample_indices=None,
    ):
        """
        Extracts read count data from the HDF5 file based on SNV/sample labels or indices.

        Parameters
        ----------
        tables : list of str, optional
            List of table types to extract. Default is ["variant", "total"].
        sources : list of str, optional
            Sample sources to filter by.
        snv_labels : list of str, optional
            Labels of SNVs to extract.
        snv_indices : list or array-like, optional
            Precomputed indices of SNVs.
        sample_labels : list of str, optional
            Labels of samples to extract.
        sample_indices : list or array-like, optional
            Precomputed indices of samples.

        Returns
        -------
        dict of {str: np.ndarray}
            A dictionary mapping table names to their extracted NumPy arrays.
        """

        # if snv_labels:

        #     snv_indices = [snv_idx[l] for l in snv_labels]

        # if sources is not None:
        #     sample_indices = self._extract_indices_by_column(
        #         "sample/metadata", "source", sources
        #     )
        # elif sample_labels is not None:

        #     sample_indices = [sample_idx[l] for l in sample_labels]

        return {
            table: self._extract_data(
                f"read_counts/{table}", snv_indices, sample_indices
            )
            for table in tables
        }

    def load_metadata(self, fname, index_name, label_col, delimiter=",", label_sep=":"):
        """
         Reads sample metadata from a CSV file, adds any new sample names to the index,
         and inserts metadata into the /sample/metadata table in the HDF5 file.

        Parameters
         ----------
        fname : str
             Path to the metadata CSV.
        index_name : str
             Name of the index to be updated.
        label_col : str | list of columns
             Column name in the CSV that contains the labels for the index.
        delimiter : str, optional
             Delimiter used in the metadata input file (default is ",").
        label_sep : str, optional
             Separator used to concatentate label column if label_col is a list (default is ":").
        """
        try:
            metadata_dict = {}
            with open(fname, newline="") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                if not isinstance(label_col, list):
                    label_col = [label_col]
                for row in reader:
                    # print(f"Row type: {type(row)}; Keys: {row.keys() if isinstance(row, dict) else row}")
                    label = label_sep.join([row[col] for col in label_col])
                    metadata_dict[label] = {
                        k: v for k, v in row.items() if k not in label_col
                    }

            if index_name not in self.global_idx:
                raise ValueError(f"Index '{index_name}' not found in global_idx.")

            self.global_idx[index_name].insert_metadata(metadata_dict)
            table_name = f"{index_name}/metadata"
            self._log_dataset_modification(table_name, "update", source_file=fname)

        except Exception:
            self.close()
            raise

    def parse_battenberg_file(self, fname, sample_label):
        """
        Parse a Battenberg file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the Battenberg file to parse.
        """
        #  Column	Description
        # chr	The chromosome of the segment
        # startpos	Start position on the chromosome
        # endpos	End position on the chromosome
        # BAF	The B-allele frequency of the segment
        # pval	P-value that is obtained when testing whether this segment should be represented by one or two states. A low p-value will result in the fitting of a second copy number state
        # LogR	The log ratio of normalised tumour coverage versus its matched normal sequencing sample
        # ntot	An internal total copy number value used to determine the priority of solutions. NOTE: This is not the total copy number of this segment!
        # nMaj1_A	The major allele copy number of state 1 from solution A
        # nMin1_A	The minor allele copy number of state 1 from solution A
        # frac1_A	Fraction of tumour cells carrying state 1 in solution A
        # nMaj2_A	The major allele copy number of state 2 from solution A. This value can be NA
        # nMin2_A	The minor allele copy number of state 2 from solution A. This value can be NA
        # frac2_A	Fraction of tumour cells carrying state 2 in solution A. This value can be NA
        # SDfrac_A	Standard deviation on the BAF of SNPs in this segment, can be used as a measure of uncertainty
        # SDfrac_A_BS	Bootstrapped standard deviation
        # frac1_A_0.025	Associated 95% confidence interval of the bootstrap measure of uncertainty
        try:
            df = pd.read_csv(fname, sep="\t")

            self.batch_sample_add([sample_label], source_file=fname)

            self.insert_sample_metadata(
                {sample_label: {"modality": Modalities.BULK.value}}
            )

            self.copy_number_bulk_view.add([sample_label])
            sample_idx = self.global_idx[GlobalIndexName.SAMPLE.value][sample_label]

            cn_labels = (
                df[["chr", "startpos", "endpos"]]
                .astype(str)
                .agg(":".join, axis=1)
                .unique()
                .tolist()
            )

            # add the copy number segments to the local index
            self.batch_copy_numbers_bulk_add(cn_labels)

            logr_dict = dict(zip(cn_labels, df["LogR"]))
            baf_dict = dict(zip(cn_labels, df["BAF"]))

            # for i, row in df.iterrows():
            #     label = cn_labels[i]

            # create the baf and logr datasets
            bafs = np.array(
                [baf_dict[segment] for segment in cn_labels], dtype="f8"
            ).reshape(-1, 1)
            logrs = np.array(
                [logr_dict[segment] for segment in cn_labels], dtype="f8"
            ).reshape(-1, 1)

            # construct the data structure for the copy number profiles
            profiles = []
            for _, row in df.iterrows():
                segments = []
                # First allele-specific state
                if (
                    pd.notna(row["nMaj1_A"])
                    and pd.notna(row["nMin1_A"])
                    and pd.notna(row["frac1_A"])
                ):
                    segments.append(
                        (
                            int(row["nMaj1_A"]),
                            int(row["nMin1_A"]),
                            float(row["frac1_A"]),
                        )
                    )
                # Optional second allele-specific state
                if (
                    pd.notna(row["nMaj2_A"])
                    and pd.notna(row["nMin2_A"])
                    and pd.notna(row["frac2_A"])
                ):
                    segments.append(
                        (
                            int(row["nMaj2_A"]),
                            int(row["nMin2_A"]),
                            float(row["frac2_A"]),
                        )
                    )
                profiles.append(np.array(segments, dtype=ALLELE_SPECIFIC_CN_DTYPE))
            profiles = np.array(
                profiles, dtype=h5py.vlen_dtype(ALLELE_SPECIFIC_CN_DTYPE)
            ).reshape(-1, 1)

            # add baf
            table_name_base = f"{SchemaGroups.COPY_NUMBERS.value}/bulk"
            seg_indices = self.indices_by_copy_numbers_bulk_label(cn_labels)
            self._add_copy_number_profiles(
                f"{table_name_base}/profile",
                profiles,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )
            self._add_copy_number_raw(
                f"{table_name_base}/baf",
                bafs,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )
            self._add_copy_number_raw(
                f"{table_name_base}/logr",
                logrs,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

        except Exception:
            self.close()
            raise

    def parse_ascat_sc_total_copy_numbers(
        self, fname, sample_label, modality=Modalities.SCDNA.value
    ):
        """
        Parse an ASCAT SC total copy number file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the ASCAT SC total copy number file to parse.
        sample_label : str
            Label for the sample being added.
        modality : str
            Modality of the sample (e.g., "scdna", "lcm"). Default is "scdna".
        """

        try:

            try:
                modality_enum = Modalities(modality)
            except ValueError:
                raise ValueError(f"Modality '{modality}' not recognized.")

            df = pd.read_csv(fname, sep="\t")

            required_columns = [
                "chromosome",
                "start",
                "end",
                "logr",
                "total_copy_number",
            ]
            missing = [col for col in required_columns if col not in df.columns]
            if missing:
                raise ValueError(f"Missing columns in {fname}: {missing}")

            self.batch_sample_add([sample_label], source_file=fname)
            self.insert_sample_metadata(
                {sample_label: {"modality": modality_enum.value}}
            )

            if modality_enum == Modalities.SCDNA:
                self.copy_number_scdna_view.add([sample_label])
            elif modality_enum == Modalities.LCM:
                self.copy_number_lcm_view.add([sample_label])

            sample_idx = self.global_idx[GlobalIndexName.SAMPLE.value][sample_label]

            seg_labels = (
                df[["chromosome", "start", "end"]]
                .astype(str)
                .agg(":".join, axis=1)
                .tolist()
            )

            # add the copy number segments to the local index
            self.batch_copy_numbers_scdna_add(seg_labels)

            logr_dict = dict(zip(seg_labels, df["logr"]))

            logrs = np.array(
                [logr_dict[segment] for segment in seg_labels], dtype="f8"
            ).reshape(-1, 1)

            profiles = np.empty(
                (df.shape[0], 1), dtype=h5py.vlen_dtype(ALLELE_SPECIFIC_CN_DTYPE)
            )
            for i, row in df.iterrows():
                total_cn = row["total_copy_number"]
                # (Ensure total_cn is valid; the check above helps here)
                profiles[i, 0] = np.array(
                    [(int(total_cn), 0, 1.0)], dtype=ALLELE_SPECIFIC_CN_DTYPE
                )

            # add baf
            table_name_base = f"{SchemaGroups.COPY_NUMBERS.value}/scdna"
            seg_indices = self.indices_by_copy_numbers_scdna_label(seg_labels)
            self._add_copy_number_profiles(
                f"{table_name_base}/profile",
                profiles,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

            self._add_copy_number_raw(
                f"{table_name_base}/logr",
                logrs,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

        except Exception:
            self.close()
            raise

    def parse_ascat_sc_allele_specific_copy_numbers(
        self, fname, sample_label, modality=Modalities.SCDNA.value
    ):
        """
        Parse an ASCAT SC total copy number file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the ASCAT SC total copy number file to parse.
        sample_label : str
            Label for the sample being added.
        modality : str
            Modality of the sample (e.g., "scdna", "lcm"). Default is "scdna".
        """

        try:

            try:
                modality_enum = Modalities(modality)
            except ValueError:
                raise ValueError(f"Modality '{modality}' not recognized.")

            df = pd.read_csv(fname, sep="\t")

            required_columns = ["chr", "startpos", "endpos", "logr", "BAF", "nA", "nB"]
            missing = [col for col in required_columns if col not in df.columns]
            if missing:
                raise ValueError(f"Missing columns in {fname}: {missing}")

            self.batch_sample_add([sample_label], source_file=fname)
            self.insert_sample_metadata(
                {sample_label: {"modality": modality_enum.value}}
            )

            if modality_enum == Modalities.SCDNA:
                self.copy_number_scdna_view.add([sample_label])
            elif modality_enum == Modalities.LCM:
                self.copy_number_lcm_view.add([sample_label])

            sample_idx = self.global_idx[GlobalIndexName.SAMPLE.value][sample_label]

            seg_labels = (
                df[["chr", "startpos", "endpos"]]
                .astype(str)
                .agg(":".join, axis=1)
                .tolist()
            )

            # add the copy number segments to the local index
            self.batch_copy_numbers_scdna_add(seg_labels)

            logr_dict = dict(zip(seg_labels, df["logr"]))

            baf_dict = dict(zip(seg_labels, df["BAF"]))

            logrs = np.array(
                [logr_dict[segment] for segment in seg_labels], dtype="f8"
            ).reshape(-1, 1)

            bafs = np.array(
                [baf_dict[segment] for segment in seg_labels], dtype="f8"
            ).reshape(-1, 1)

            profiles = np.empty(
                (df.shape[0], 1), dtype=h5py.vlen_dtype(ALLELE_SPECIFIC_CN_DTYPE)
            )
            for i, row in df.iterrows():
                nA = int(row["nA"])
                nB = int(row["nB"])
                total_cn = row["total_copy_number"]
                # (Ensure total_cn is valid; the check above helps here)
                profiles[i, 0] = np.array(
                    [(nA, nB, 1.0)], dtype=ALLELE_SPECIFIC_CN_DTYPE
                )

            # add baf
            table_name_base = f"{SchemaGroups.COPY_NUMBERS.value}/scdna"
            seg_indices = self.indices_by_copy_numbers_scdna_label(seg_labels)
            self._add_copy_number_profiles(
                f"{table_name_base}/profile",
                profiles,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

            self._add_copy_number_raw(
                f"{table_name_base}/logr",
                logrs,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

            self._add_copy_number_raw(
                f"{table_name_base}/baf",
                bafs,
                seg_indices,
                [sample_idx],
                source_file=fname,
            )

        except Exception:
            self.close()
            raise

    def segment_lookup(self, locus, group_name):
        """
        Given a genomic locus, return the segment index for the specified group.

        Parameters
        ----------
        locus : str
            Genomic locus in the format "chr:pos".
        group_name : str
            Name of the group to search for the locus.

        Returns
        -------
        int
            The segment index corresponding to the locus.
        """
        chrom, start = locus.split(":")
        if group_name == Modalities.BULK.value:
            index = self.get_copy_numbers_bulk_index()
        elif group_name == Modalities.LCM.value:
            index = self.get_copy_numbers_lcm_index()
        elif group_name == Modalities.SCDNA.value:
            index = self.get_copy_numbers_scdna_index()

        for label, idx in index.items():
            seg_chrom, seg_start, seg_end = label.split(":")
            if seg_chrom == chrom:
                if int(start) >= int(seg_start) and int(start) <= int(seg_end):
                    return label, idx
        return None, None

    def _add_copy_number_raw(
        self, table_name, dat, seg_indices, sample_indices, source_file=""
    ):
        """
        Add 2-d copy number data (logr, baf) to the HDF5 dataset.
        """

        sample_indices = np.asarray(sample_indices)

        # Ensure sorting

        sample_sort = np.argsort(sample_indices)

        # seg_indices_sorted = seg_indices[seg_sort]
        sample_indices_sorted = sample_indices[sample_sort]

        # # Sort data to match HDF5 access pattern
        dat_sorted = dat[:, sample_sort]

        # Ensure dat is a NumPy array and validate shape

        if dat.shape != (len(seg_indices), len(sample_indices)):
            raise ValueError(
                f"Shape of `dat` {dat.shape} does not match seg_indices {len(seg_indices)} x sample_indices {len(sample_indices)}"
            )

        # Safe assignment with index alignment
        for i, seg in enumerate(seg_indices):
            self[table_name][seg, sample_indices_sorted] = dat_sorted[i, :]

        self._log_dataset_modification(
            table_name, operation="update", source_file=source_file
        )

    def _add_copy_number_profiles(
        self, table_name, dat, seg_indices, sample_indices, source_file=""
    ):
        """
        Addy  to the HDF5 dataset.
        """

        if dat.shape != (len(seg_indices), len(sample_indices)):
            raise ValueError(
                f"Shape of `dat` {dat.shape} does not match seg_indices {len(seg_indices)} x sample_indices {len(sample_indices)}"
            )

        for i, seg in enumerate(seg_indices):
            for j, sample in enumerate(sample_indices):

                self[table_name][seg, sample] = np.asarray(
                    dat[i, j], dtype=ALLELE_SPECIFIC_CN_DTYPE
                )

        self._log_dataset_modification(
            table_name, operation="update", source_file=source_file
        )

    def parse_pyclone_file(self, source_file):
        """
        Parse a PyClone standard output file and extract the SNV cluster assignments.
        These will update the SNV global index cluster assignments.

        Parameters
        ----------
        source_file : str
            Path to the PyClone file to parse.



        # PYCLONE OUTPUT FORMAT:
        # The results file output by write-results-file is in tab delimited format. There six columns:

        # mutation_id - Mutation identifier as used in the input file.

        # sample_id - Unique identifier for the sample as used in the input file.

        # cluster_id - Most probable cluster or clone the mutation was assigned to.

        # cellular_prevalence - Proportion of malignant cells with the mutation in the sample.
        # This is also called cancer cell fraction (CCF) in the literature.

        # cellular_prevalence_std - Standard error of the cellular_prevalence estimate.

        # cluster_assignment_prob - Posterior probability the mutation is assigned to the cluster.
        # This can be used as a confidence score to remove mutations with low probability of belonging to a cluster.

        """

        # TODO: What do do about CCFs and cluster assignment probs?
        pyclone = pd.read_table(source_file, sep="\t")

        if self.verbose:
            print(
                f"#Parsing PyClone output file {source_file} with {pyclone.shape[0]} rows."
            )

        # Extract the SNV labels
        snv_labels = pyclone["mutation_id"].tolist()
        sample_labels = pyclone["sample_id"].tolist()

        _ = self.batch_snv_add(snv_labels, source_file=source_file)

        _ = self.batch_sample_add(sample_labels, source_file=source_file)

        cluster_dict = dict(zip(snv_labels, pyclone["cluster_id"].astype(int)))

        # update the clusters
        self.update_snv_clusters(cluster_dict, source_file=source_file)

        if self.verbose:
            print(f"#{len(cluster_dict)} SNV cluster assignments updated.")

        # self._log_dataset_modification("index/SNV/cluster", operation="update", source_file=source_file)

    def initialize_pseudobulk_layer(self, sample_label, source_file=""):
        """
        Automatically generates copynumber datasets for a pseudob-bulk analysis and
        adds sample to the sample index.
        """
        pseudo_bulk_schema = {sample_label: COPY_NUMBER_LAYER_DICT}
        self.add_sample(sample_label, data={"source": "pseudo-bulk"}, overwrite=True)
        self._recursive_build(pseudo_bulk_schema, "copy_numbers")
        self._log_dataset_modification(
            f"copy_numbers/{sample_label}", operation="create", source_file=source_file
        )

    def _is_connected(self):
        if not hasattr(self, "file"):
            return False
        else:
            return self.file.id.valid
        # self.file = h5py.File(self.filename, "a")  # Append mode

    def close(self):
        """
        Close the HDF5 file safely.

        This should be called when finished working with the DNAStream object.
        """
        for _, global_idx in self.global_idx.items():
            global_idx.save_index()

        if self._is_connected():

            self.file.close()
            if self.verbose:
                print(f"#Stream to connection {self.filename} closed.")
        else:
            print("#Stream to connection is already closed.")
