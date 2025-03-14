import sys
import os
import getpass
import socket
import pathlib
import time
from datetime import datetime
import functools
import pandas as pd
import h5py
import numpy as np

# import json
from .index_manager import BaseIndex, GlobalIndex

from .schema import (
    SCHEMA,
    STRUCT_ARRAYS,
    META_TABLES,
    MODALITIES,
    COPY_NUMBER_LAYER_DICT,
)
from .datatypes import EDGE_LIST_DTYPE


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

    BULK = "bulk"
    LCM = "lcm"
    SCDNA = "scdna"

    SNV = "SNV"
    SNP = "SNP"
    SAMPLE = "sample"
    TREE = "tree"
    CNA = "CNA"
    CLONAL = "clonal"

    TREE_TYPES = ["SNV", "CNA"]

    GLOBAL_INDICES = [SNV, SAMPLE, SNP]
    LOCAL_INDICES = [f"tree/{ttype}_trees" for ttype in TREE_TYPES] + [
        f"copy_numbers/{mod}" for mod in MODALITIES
    ]
    TREES = [SNV, CNA, CLONAL]

    def __init__(
        self, filename, initialize=True, verbose=False, id=None, sex=None, safe=True
    ):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose
        self.safe = safe

        self.file = h5py.File(filename, "a")  # Append mode (does not overwrite)
        if self.verbose:
            print(f"#Stream to connection {self.filename} open...")

        # only for 1D unchunked data, like logs, dataframes, tree lists
        # multi-dimensional data tables like READ_COUNTS, must be built separately
        # to optimize chunking and shape specification
        if initialize:
            self._recursive_build(SCHEMA)

        self.global_idx = {
            name: GlobalIndex(
                self.file,
                f"index/{name}",
                metadata_dtype=SCHEMA["index"][name]["metadata"]["dtype"],
            )
            for name in SCHEMA["index"].keys()
        }
        self.local_idx = {
            name: BaseIndex(self.file, name, metadata_dtype=self._local_idx_dtype(name))
            for name in DNAStream.LOCAL_INDICES
        }
        # Create group structure for read counts
        # for key, _ in READ_COUNTS.items():
        #     # self.add_dataset_to_file(
        #     #     key,
        #     #     shape=(0, 0),
        #     #     maxshape=(None, None),
        #     #     dtype=dtype,
        #     #     compression="gzip",
        #     #     chunks=(1, 5000),
        #     # )
        # TODO: add dim labels
        #     self.file[key].dims[0].label = DNAStream.SNV
        #     self.file[key].dims[1].label = DNAStream.SAMPLE

        if id:
            self.set_patient_id(id)
        else:
            self.set_patient_id("")

        if sex:
            self.set_patient_sex(sex)
        else:
            self.set_patient_sex("")

        if not self.safe and id:
            self.set_patient_id(id)
        if not self.safe and sex:
            self.set_patient_sex(sex)

    def __str__(self):
        """To string method"""
        m = self.file[f"{DNAStream.SNV}/label"].shape[0]
        n = self.file[f"{DNAStream.SAMPLE}/label"].shape[0]

        mystr = f"DNAStream object with {m} SNVs and {n} samples"
        mystr += f"\nHDF5 File: {self.filename}"
        mystr += f"\nPatient: {self.file.attrs['id']}, sex: {self.file.attrs['sex']}"

        return mystr

    def _dtype(self, table):
        return self.file[table].dtype

    def _local_idx_dtype(self, table):
        if f"{table}/metadata" not in self.file:
            print(f"{table}/metadata")
        dtype = self.file[f"{table}/metadata"].dtype
        return dtype

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
                new_path = f"{path}/{key}" if path else key  # Construct the HDF5 path

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
                print(f"Creating dataset {path}...")
            self.file.create_dataset(path, dtype=dtype, **kwargs)
            if columns:
                self.file[path].attrs["columns"] = columns
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
        log = self.file["metadata/log"]
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

    @timeit
    def add_read_counts(self, fname, source, location=None):
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
            rc = pd.read_csv(
                fname, names=["snv", "sample", "var", "total"], header=None, skiprows=1
            )

            samples = rc["sample"].unique().tolist()

            snvs = rc["snv"].unique().tolist()

            snv_idx = self.batch_add_snvs(snvs, source_file=fname)

            sample_idx = self.batch_add_samples(samples, source_file=fname)
            sample_indices = list(sample_idx.values())
            if source:
                self._update_value(
                    sample_indices, f"{DNAStream.SAMPLE}/metadata", "source", source
                )

            if location:
                self._update_value(
                    sample_indices, f"{DNAStream.SAMPLE}/metadata", "location", location
                )

            # Map indices for SNVs and samples
            rc["snv_idx"] = rc["snv"].map(snv_idx)
            rc["sample_idx"] = rc["sample"].map(sample_idx)

            # Convert to NumPy arrays for efficient updates
            snv_indices_arr = rc["snv_idx"].to_numpy()
            sample_indices_arr = rc["sample_idx"].to_numpy()
            var_counts = rc["var"].to_numpy()
            total_counts = rc["total"].to_numpy()

            # **Sort indices to satisfy HDF5 fancy indexing rules**
            sorted_order = np.lexsort((sample_indices_arr, snv_indices_arr))
            snv_indices_arr = snv_indices_arr[sorted_order]
            sample_indices_arr = sample_indices_arr[sorted_order]
            var_counts = var_counts[sorted_order]
            total_counts = total_counts[sorted_order]

            # Update dataset in batch instead of looping
            dat = self.file[f"read_counts"]
            unique_snv_indices = np.unique(snv_indices_arr)  # Get unique SNV indices

            for snv in unique_snv_indices:
                mask = snv_indices_arr == snv  # Select all entries for this snv index
                dat["variant"][snv, sample_indices_arr[mask]] = var_counts[mask]
                dat["total"][snv, sample_indices_arr[mask]] = total_counts[mask]

            for arr in ["variant", "total"]:
                self._log_dataset_modification(
                    f"read_counts/{arr}", operation="update", source_file=fname
                )

        except Exception as e:
            self.close()
            raise Exception(e)

    def create_global_index(self, name, metadata_dtype):
        """
        Add an empty global index to the HDF5 file.

        Parameters
        ----------
        name : str
            Name of the index.
        metadata_dtype : numpy.dtype
            Data type of the index metadata.
        """
        self.global_idx[name] = GlobalIndex(self.file, name, metadata_dtype)

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

        dat = self.file[f"read_counts"]
        dat["variant"][snv_idx, sample_idx] = var
        dat["total"][snv_idx, sample_idx] = total

    def snv_label_to_idx(self, label):
        """
        Retrieve the index of an SNV label.

        Parameters
        ----------
        label : str
            SNV label to look up.

        Returns
        -------
        int or None
            Index of the SNV if found, else None.
        """
        return self.global_idx[DNAStream.SNV].sample_label_to_idx(label)

    def sample_label_to_idx(self, label):
        """
        Retrieve the index of a sample label.

        Parameters
        ----------
        label : str
            sample label to look up.

        Returns
        -------
        int or None
            Index of the sample if found, else None.
        """
        return self.global_idx[DNAStream.SAMPLE][label]

    def _resize_all(self, m=None, n=None):
        """
        Resize the HDF5 datasets to accommodate additional SNVs or samples.

        Parameters
        ----------
        m : int, optional
            New size for the SNV index (default is None, which does not resize).
        n : int, optional
            New size for sample index (default is None, which does not resize).
        """

        if m:
            for snv_data in META_TABLES:
                if snv_data == "log":
                    continue
                group_path = f"{DNAStream.SNV}/{snv_data}"

                self.file[group_path].resize((m,))
        else:
            m = self.global_idx["SNV"].size()

        if n:
            for sample_data in META_TABLES:
                if sample_data == "log":
                    continue
                group_path = f"{DNAStream.SAMPLE}/{sample_data}"
                self.file[group_path].resize((n,))
        else:
            n = self.global_idx["sample"].size()

        group_path = "read_counts"
        for reads in ["variant", "total"]:

            mat = self.file[f"{group_path}/{reads}"]
            mat.resize((m, n))

    def add_snv(self, label, cluster=None, data=None, overwrite=False):
        """
        Add a single SNV (Single Nucleotide Variant) to the dataset.

        Parameters
        ----------
        label : str
            The SNV label (e.g., concatenated chromosome, position, reference, and alternate allele).
        cluster : int, optional
            Cluster ID associated with the SNV (default is None).
        data : dict, optional
            Additional metadata to store for the SNV (default is None).
        overwrite : bool, optional
            If True, overwrite existing SNV data if it already exists (default is False).

        Returns
        -------
        int
            The assigned index of the SNV in the dataset.
        """
        return self.add_item(
            label, index=DNAStream.SNV, cluster=cluster, data=data, overwrite=overwrite
        )

    def add_sample(self, label, cluster=None, data=None, overwrite=False):
        """
        Add a single sample to the dataset.

        Parameters
        ----------
        label : str
            The sample label (e.g., patient or sample ID).
        cluster : int, optional
            Cluster ID associated with the sample (default is None).
        data : dict, optional
            Additional metadata to store for the sample (default is None).
        overwrite : bool, optional
            If True, overwrite existing sample data if it already exists (default is False).

        Returns
        -------
        int
            The assigned index of the sample in the dataset.
        """
        return self._add_item(
            label,
            index=DNAStream.SAMPLE,
            cluster=cluster,
            data=data,
            overwrite=overwrite,
        )

    def _add_item(
        self,
        label: str,
        index_name,
        index_dict=None,
        cluster=None,
        data=None,
        overwrite=False,
    ):
        """
        Internal method to add an item (SNV or sample) to the dataset.

        Parameters
        ----------
        label : str
            The unique identifier for the SNV or sample.
        index_name : str
            The index type, either DNAStream.SNV or DNAStream.SAMPLE.
        index_dict : dict, optional
            Pre-loaded index dictionary for quick lookups (default is None).
        cluster : int, optional
            Cluster ID for the SNV or sample (default is None).
        data : dict, optional
            Additional metadata to associate with the SNV or sample (default is None).
        overwrite : bool, optional
            If True, overwrite existing entry (default is False).

        Returns
        -------
        int
            The assigned index of the item in the dataset.

        Raises
        ------
        ValueError
            If the label already exists and `overwrite=False`, a warning is printed, and no action is taken.
        """
        label = label.encode("utf-8")
        idx = self.idx_by_label(label, index_name, index_dict)

        if not idx:
            new_idx = self.file[f"{index_name}/label"].shape[0]
            if index_name == DNAStream.SNV:
                self._resize_all(m=new_idx + 1)
            else:
                self._resize_all(n=new_idx + 1)
        else:
            if overwrite:
                new_idx = idx
            else:
                print(
                    f"{label} exists in {index_name}, use overwrite=True to overwrite metadata."
                )
                return idx

        self.file[f"{index_name}/label"][new_idx] = label
        self.file[f"{index_name}/index"][label] = new_idx
        if data:
            self.file[f"{index_name}/metadata"][new_idx] = data
        if cluster:
            self.file[f"{index_name}/cluster"][new_idx] = cluster

        return new_idx

    def load_snv_index(self):
        """
        Load the SNV index from the HDF5 file.

        Returns
        -------
        dict
            A dictionary mapping SNV labels to their respective indices.
        """
        return self.global_idx[DNAStream.SNV].get_index()

    def load_sample_index(self):
        """
        Load the sample index from the HDF5 file.

        Returns
        -------
        dict
            A dictionary mapping sample labels to their respective indices.
        """
        return self.global_idx[DNAStream.SNV].get_index()

    def batch_add_snvs(self, labels, source_file=""):
        """
        Batch add multiple SNVs to the dataset.

        Parameters
        ----------
        labels : list of str
            List of SNV labels to add.
        source_file : str, optional
            Path to the source file from which the SNVs are being added (default is an empty string).

        Returns
        -------
        dict
            A dictionary mapping SNV labels to their respective indices.
        """
        snv_dict = self.global_idx[DNAStream.SNV].add(labels, source_file)

        self._resize_all(m=self.global_idx[DNAStream.SNV].size())

        return snv_dict

    def batch_add_samples(self, labels, source_file="", metadata=None):
        """
        Batch add multiple samples to the dataset.

        Parameters
        ----------
        labels : list of str
            List of sample labels to add.
        source_file : str, optional
            Path to the source file from which the samples are being added (default is an empty string).

        Returns
        -------
        dict
            A dictionary mapping sample labels to their respective indices.
        """
        sample_dict = self.global_idx[DNAStream.sample].add(labels, source_file)

        self._resize_all(n=self.global_idx[DNAStream.sample].size())

        return sample_dict

    def get_snv_data(self, indices=None):
        """
        Retrieve SNV data from the HDF5 file.

        Parameters
        ----------
        indices : list of int, optional
            A list of SNV indices to retrieve. If None, retrieves all SNVs.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the SNV data.
        """
        return self._get_data(dataset_name=f"{DNAStream.SNV}/metadata", indices=indices)

    def get_snv_log(self):
        """
        Retrieve the SNV index update log.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for SNV updates.
        """
        return self._get_data(dataset_name=self.global_idx[DNAStream.SNV].log_name)

    def get_sample_log(self):
        """
        Retrieve the sample index update log.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for sample updates.
        """
        return self._get_data(dataset_name=self.global_idx[DNAStream.sample].log_name)

    def get_dataset_log(self):
        """
        Retrieve the dataset log

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for SNV updates.
        """
        return self._get_data(dataset_name=f"metadata/log")

    def get_sample_data(self, indices=None):
        """
        Retrieve sample data from the HDF5 file.

        Parameters
        ----------
        indices : list of int, optional
            A list of sample indices to retrieve. If None, retrieves all samples.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the sample data.
        """
        return self._get_data(
            dataset_name=f"{DNAStream.SAMPLE}/metadata", indices=indices
        )

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
        dataset = self.file[dataset_name]

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
            dataset_name=f"{DNAStream.SNV}/metadata",
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
            dataset_name=f"{DNAStream.SAMPLE}/metadata",
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
        dataset = self.file[dataset_name]
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
        dataset = self.file[dataset_name]
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

            snv_idx = self.batch_add_snvs(snv_labels, source_file=fname)

            # sort dataframe according to the newly assigned indices in DNAStream
            maf.loc[:, "snv_idx"] = maf["label"].map(snv_idx)
            maf = maf.sort_values("snv_idx")
            indices = maf["snv_idx"].tolist()
            snv_data = maf[[val for _, val in column_dict.items()]]

            self.add_snv_data(indices, snv_data, source_file=fname)

        except Exception as e:

            self.close()
            raise Exception(e)

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
        for table in table_list:
            old_size = self.file[table].shape[0]
            new_size = old_size + n
            self.file[table].resize((new_size,))

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
        except Exception as e:
            self.close()
            raise Exception(e)

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
        table_name = f"{DNAStream.TREE}/{tree_type}_trees"
        tree_list = wrap_list(tree_list)

        tree_dict = self.local_idx[table_name].resize(len(tree_list), prefix="snv_tree")

        tree_indices = list(tree_dict.values())

        numtrees = len(tree_dict)
        self._expand([f"{table_name}/trees"], numtrees)

        # Add the tree metadata to the file
        data_dtype = self._dtype(f"{DNAStream.TREE}/{tree_type}_trees/metadata")
        tree_dtype = self._dtype(f"{DNAStream.TREE}/{tree_type}_trees/trees")

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

        self.file[f"{table_name}/metadata"][tree_indices] = dat

        tree_dtype = SCHEMA[DNAStream.TREE][f"{tree_type}_trees"]["trees"]["dtype"]
        tree_structured = np.array(
            [np.array(tree, dtype=EDGE_LIST_DTYPE) for tree in tree_list],
            dtype=tree_dtype,
        )
        self.file[f"{table_name}/trees"][tree_indices] = tree_structured

        self._log_dataset_modification(
            f"{DNAStream.TREE}/{tree_type}_trees", "update", source_file=source_file
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

            dataset_name = f"{DNAStream.TREE}/{tree_type}_trees/metadata"

            if self.safe and not self._check_safe(dataset_name, source_file):
                print(
                    "#Warning! Attempting overwrite in Safe mode, use safe=F, to force append trees."
                )
                return

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

        except Exception as e:
            self.close()
            raise Exception(e)

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
        data_dtype = self._dtype(f"{DNAStream.TREE}/SNV_trees/metadata")

        dat = np.array([("", method, score, rank, "")], dtype=data_dtype)
        self._add_trees(edge_list, "SNV", data=dat)
        self._log_dataset_modification(f"{DNAStream.TREE}/SNV_tree", operation="update")

    def add_pyclone_file(self, fname):
        pass

    def _extract_indices_by_column(self, dataset_name, name, values):
        vals = self.file[dataset_name][name][:]  # Load the column data
        indices = np.where(np.isin(vals, values))[0]  # Get matching indices
        return indices

    def load_indices(self):
        """
        Wrapper to return both SNV and sample indices as dictionaries
        """
        return self.load_snv_index(), self.load_sample_index()

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
            arr = self.file[dataset_name][snv_indices, :]
            return arr[:, sample_indices]
        elif snv_indices:
            return self.file[dataset_name][snv_indices]
        elif sample_indices:
            return self.file[dataset_name][:, sample_indices]
        else:
            return self.file[dataset_name][:]

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

        snv_idx, sample_idx = self.load_indices()
        if snv_labels:

            snv_indices = [snv_idx[l] for l in snv_labels]

        if sources is not None:
            sample_indices = self._extract_indices_by_column(
                "sample/metadata", "source", sources
            )
        elif sample_labels is not None:

            sample_indices = [sample_idx[l] for l in sample_labels]

        return {
            table: self._extract_data(
                f"read_counts/{table}", snv_indices, sample_indices
            )
            for table in tables
        }

    def add_copy_numbers(self, labels, values, source_file=""):
        """
        Add copy number data to the HDF5 dataset.

        Parameters
        ----------
        labels : list of str
            List of copy number labels.
        values : list of float
            List of copy number values.
        source_file : str, optional
            Path to the source file from which the copy numbers are being added (default is an empty string).

        Raises
        ------
        Exception
            If there is an error processing the file or adding copy numbers.

        Notes
        -----
        - The function checks if copy numbers from the source file already exist.
        - Logs dataset modifications after adding copy numbers.
        """
        try:
            source_file = str(pathlib.Path(source_file).resolve())

            dataset_name = f"{DNAStream.COPY_NUMBER}/metadata"

            if self.safe and not self._check_safe(dataset_name, source_file):
                print(
                    "#Warning! Attempting overwrite in Safe mode, use safe=F, to force append copy numbers."
                )
                return

            copy_number_dict = self.local_idx[dataset_name].resize(
                len(labels), prefix="copy_number"
            )

            copy_number_indices = list(copy_number_dict.values())

            num_copy_numbers = len(copy_number_dict)
            self._expand([f"{dataset_name}/metadata"], num_copy_numbers)

            # Add the copy number metadata to the file
            data_dtype = self._dtype(f"{DNAStream.COPY_NUMBER}/metadata")
            # SCHEMA[DNAStream.COPY_NUMBER]["data"]["dtype"]
            dat = np.array(
                [(lab, val, source_file) for lab, val in zip(copy_number_dict, values)],
                dtype=data_dtype,
            )

            self.file[f"{dataset_name}"][copy_number_indices] = dat

            self._log_dataset_modification(
                f"{DNAStream.COPY_NUMBER}/metadata", "update", source_file=source_file
            )

        except Exception as e:
            self.close()
            raise Exception(e)

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

            self.add_sample(sample_label, data={"source": "battenberg"}, overwrite=True)

            # cn_labels = df[['chr', 'startpos', 'endpos']].astype(str).agg(':'.join, axis=1).unique().tolist()

            # self.lo

            # self.add_copy_numbers(cn_labels, cn_values, source_file=fname)

        except Exception as e:
            self.close()
            raise Exception(e)

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

    def close(self):
        """
        Close the HDF5 file safely.

        This should be called when finished working with the DNAStream object.
        """
        for _, global_idx in self.global_idx.items():
            global_idx.save_index()

        self.file.close()
        if self.verbose:
            print(f"#Stream to connection {self.filename} closed.")
