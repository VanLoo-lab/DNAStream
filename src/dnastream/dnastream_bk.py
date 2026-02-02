import getpass
import socket
import pathlib

from datetime import datetime
import functools
import pandas as pd
import h5py
import numpy as np
from enum import Enum

# impot mixins
from .io import IOMixin

from .utils import (
    wrap_list,
    # full_path,
    # timeit
)

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


"""
DNAStream is opinionated about biology, but boringly consistent about access.
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


class DNAStream(IOMixin):
    """
    DNAStream is an HDF5-based data structure for efficient storage, indexing,
    and retrieval of processed DNA sequencing data across multiple sequencing modalities.

    It provides structured storage for SNV read counts, copy number profiles, and metadata
    from bulk, LCM, and single-cell sequencing. The design ensures consistency across
    different data modalities and enables efficient querying and updating.



    """

    def __init__(self, filename, verbose=False, id=None, sex=None, safe=True):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose

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
                "Warning: You are connecting outside of a context manager. This is not recommended as runtime errors may corrupt the file."
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
            idx_view_mapping_table="copy_numbers/lcm/idx_view_mapping",
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
            idx_view_mapping_table="copy_numbers/bulk/idx_view_mapping",
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
            idx_view_mapping_table="copy_numbers/scdna/idx_view_mapping",
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
        dtype,
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
