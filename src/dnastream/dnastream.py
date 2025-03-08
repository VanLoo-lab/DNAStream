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
import json 

from .schema import SCHEMA, STRUCT_ARRAYS, META_TABLES, MODALITIES, READ_COUNTS
from .datatypes import EDGE_LIST_DTYPE 
#TODO
# - add pyclone addition
# - add patient metadata
# - add bulk phylogenies as edge lists

"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                #dataframe structure containing quality scores, number of callers, etc
     |── cluster
     |── index_map             #json string for fast loading and saving
     |-- log
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
     |-- log
 |-- trees/
 |   |-- SNV_trees/  
 |   |     |- trees      #edge lists of clusters
 |   |     |- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |- index_map  #hold the index map
 |   |-- CNA_trees     
 |   |     |- trees      #edge lists (*) probably changed to Newick strings
 |   |     |- data      #holds the label, likelihood, rank, method used to generate, etc
 |   |     |- index_map
 |   |-- clonal_trees (joint CNA SNV tree)/
 |   |     |-- tree 
 |   |     |-- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |-- index_map
 |   |     |-- genotypes (*) (structured array of node/snv/x/y/x_bar/y_bar)
 |   |     |-- clonal proportions (*) (U)
 |   |     |-- sample assignment (*)
 |   |    
 |-- copy_number/ (*)
 |   ├── /bulk/
 |   │   ├── /segments    # Bulk-specific segment index
 |   │   ├── /profiles       # Tensor: (sample, segment, allele CN, proportion μ)
 |   │   ├── /metadata    # Bulk-specific metadata
 |   |   |-- /log          # logging
 |   │
 |   ├── /single_cell/
 |   │   ├── /segments    # Single-cell segment index
 |   │   ├── /profiles       # (sample, segment) → allele CN tuple
 |   │   ├── /metadata    # Single-cell metadata
 |   │
 |   ├── /lcm/
 |   │   ├── /segments    # LCM-specific segment index
 |   │   ├── /profiles       # (sample, segment) → allele CN tuple
 |   |   |--/logR
 |   |   |--/baf
 |   │   ├── /metadata    # LCM-specific metadata
 ├── metadata/ 
     |-- log                  # Metadata storage
 │   ├── sample_info  (*)            # Sample IDs
 │   ├── processing_parameters (*)
"""




def timeit(func):
    """Decorator to measure execution time of a function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()  # Start timer
        result = func(*args, **kwargs)    # Run the function
        end_time = time.perf_counter()    # End timer
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
    SAMPLE = "sample"
    TREE = "tree"
    CNA = "CNA"
    CLONAL= "clonal"
    

    INDICES = [SNV, SAMPLE]
    TREES = [SNV, CNA, CLONAL ]

    def __init__(self, filename, initialize=True, verbose=False):

        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose
    
    
        self.file = h5py.File(filename, "a")  # Append mode (does not overwrite)
        if self.verbose:
            print(f"#Stream to connection {self.filename} open...")
        
        
        #only for 1D unchunked data, like logs, dataframes, tree lists
        #multi-dimensional data tables like READ_COUNTS, must be built separately 
        # to optimize chunking and shape specification
        if initialize:
            self._recursive_build(SCHEMA)



            # Create group structure for read counts
            for key, dtype in READ_COUNTS.items():
                        self.add_dataset_to_file(key, shape=(0,0), maxshape=(None,None), dtype=dtype, 
                                                    compression="gzip", chunks=(1, 5000))
                        self.file[key].dims[0].label = DNAStream.SNV
                        self.file[key].dims[1].label = DNAStream.SAMPLE                      


        
    def __str__(self):
        """To string method"""
        m = self.file[f"{DNAStream.SNV}/label"].shape[0]
        n = self.file[f"{DNAStream.SAMPLE}/label"].shape[0]

        mystr = f"DNAStream object with {m} SNVs and {n} samples" 
        mystr += f"\nHDF5 File: {self.filename}"
      
        return mystr
    

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
            for key, datasets in schema.items():
                new_path = f"{path}/{key}" if path else key  # Handle root case
                
                if isinstance(datasets, dict):  # If it's another dictionary, recurse
                    self._recursive_build(datasets, new_path)
                else:  # If it's a dataset, create it
                    dtype = datasets  # Since `datasets` holds dtype here
                    
                    columns = []
                    if key in STRUCT_ARRAYS:
                        columns = list(dtype.names)  # Get column names from dtype
                    
                 
                    self.add_dataset_to_file(
                        new_path, shape=(0,), maxshape=(None,),
                        dtype=dtype, compression="gzip",
                        columns=columns
                    )

    

    def add_dataset_to_file(self, path,  
                             dtype=h5py.string_dtype("utf-8"), 
                             columns=[], 
                             source_file="",
                             **kwargs):
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
                self.file[path].attrs['columns'] = columns
            self._log_dataset_modification(path, operation="create", source_file=source_file)
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
            timestamp_str, dataset_name.encode("utf-8"), operation.encode("utf-8"),
            user, hostname, source_file
        )

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
            rc = pd.read_csv(fname, names=["snv","sample", "var", "total"], header=None, skiprows=1)
       

            samples = rc["sample"].unique().tolist()

            snvs = rc["snv"].unique().tolist()

            snv_idx = self.batch_add_snvs(snvs, source_file=fname)

            sample_idx = self.batch_add_samples(samples, source_file=fname)
            sample_indices = list(sample_idx.values())
            if source:
                self._update_value(sample_indices, f"{DNAStream.SAMPLE}/data", "source", source)
            
            if location:
                self._update_value(sample_indices, f"{DNAStream.SAMPLE}/data", "location", location)

            

                
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
            dat = self.file[f"read_counts/{source}"]
            unique_snv_indices = np.unique(snv_indices_arr)  # Get unique SNV indices
            
    

            for snv in unique_snv_indices:
                mask = snv_indices_arr == snv  # Select all entries for this snv index
                dat["variant"][snv, sample_indices_arr[mask]] = var_counts[mask]
                dat["total"][snv, sample_indices_arr[mask]] = total_counts[mask]
            
            for arr in ["variant", "total"]:
                self._log_dataset_modification(f"read_counts/{source}/{arr}", operation="update", source_file=fname)





        except Exception as e:
            self.close()
            raise Exception(e)





        
    def _add_read_count(self, source,snv_idx, sample_idx, var=0, total=0):
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
    
        dat = self.file[f"read_counts/{source}"]
        dat['variant'][snv_idx,sample_idx] = var 
        dat['total'][snv_idx,sample_idx] = total
    

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
        return self._idx_by_label(label, index_name=DNAStream.SNV)
    
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
        return self._idx_by_label(label, index_name=DNAStream.SAMPLE)
    
    def _idx_by_label(self, label, index_name):
        """
        Internal function to retrieve the index from a label.

        Parameters
        ----------
        label : str
            Label to look up.
        index_name : str
            Name of the index to use.


        Returns
        -------
        int or None
            Index if found, else None.
        """
        labs = self.file[f"{index_name}/label"]

        if len(labs) ==0:
            return None

        indices = np.where(labs[:] ==label)[0]
        if len(indices) ==0:
            return None 
        elif len(indices) == 1:
            return indices[0]
        else:
            self.close()
            raise ValueError("label is associated with multiple indices and could not be added.")
     
    


    def _label_by_idx(self, idx, index_name):
        """
        Internal function to label the index from an index.

        Parameters
        ----------
        idx : int
            index to look up.
        index_name : str
            Name of the index to use.


        Returns
        -------
        str 
            label if found
        """
        try:
            label = self.file[f"{index_name}/label"][idx] 
        except Exception as e:
            self.close()
            raise Exception(e)

        return label   




    def _resize_all(self, m=None,n=None ):
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
            m= self.file[f"{DNAStream.SNV}/label"].shape[0]

        if n:
            for sample_data in META_TABLES:
                if sample_data == "log":
                    continue
                group_path = f"{DNAStream.SAMPLE}/{sample_data}"
                self.file[group_path].resize((n,))
        else:
            n = self.file[f"{DNAStream.SAMPLE}/label"].shape[0]


        for modality in MODALITIES:
            group_path = f"read_counts/{modality}"
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
        return self.add_item(label, index=DNAStream.SNV, cluster=cluster, data=data, overwrite=overwrite)


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
        return self.add_item(label, index=DNAStream.SAMPLE, cluster=cluster, data=data, overwrite=overwrite)


    def _add_item(self, label: str, index_name, index_dict=None, cluster=None, data=None, overwrite=False):
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
                print(f"{label} exists in {index_name}, use overwrite=True to overwrite metadata.")
                return idx

        self.file[f"{index_name}/label"][new_idx] = label
        self.file[f"{index_name}/index_map"][label] = new_idx
        if data:
            self.file[f"{index_name}/data"][new_idx] = data
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
        return self._load_index(DNAStream.SNV)


    def load_sample_index(self):
        """
        Load the sample index from the HDF5 file.

        Returns
        -------
        dict
            A dictionary mapping sample labels to their respective indices.
        """
        return self._load_index(DNAStream.SAMPLE)


    def save_snv_index(self, index_dict):
        """
        Save the SNV index to the HDF5 file.

        Parameters
        ----------
        index_dict : dict
            Dictionary mapping SNV labels to their respective indices.
        """
        return self._save_index(index_dict, f"{DNAStream.SNV}/index_map")


    def save_sample_index(self, index_dict):
        """
        Save the sample index to the HDF5 file.

        Parameters
        ----------
        index_dict : dict
            Dictionary mapping sample labels to their respective indices.
        """
        return self._save_index(index_dict, f"{DNAStream.SAMPLE}/index_map")


    def _load_index(self, index_name):
        """
        Load the label-to-index mapping from the HDF5 file.

        Parameters
        ----------
        index_name : str
            The index type, either DNAStream.SNV or DNAStream.SAMPLE.

        Returns
        -------
        dict
            A dictionary mapping labels to indices.
        """
        if index_name in self.file:
            index_data = self.file[index_name][()]
            return json.loads(index_data[0]) if len(index_data) > 0 else {}
        return {}


    def _save_index(self, index_dict, index_name):
        """
        Save the label-to-index mapping into the HDF5 file.

        Parameters
        ----------
        index_dict : dict
            Dictionary mapping labels to indices.
        index_name : str
            The index type, either DNAStream.SNV or DNAStream.SAMPLE.
        """
        index_json = json.dumps(index_dict)  # Convert dictionary to JSON string
        self.file[index_name].resize((1,))  # Ensure space in dataset
        self.file[index_name][0] = index_json  # Store JSON string in HDF5


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
        return self._batch_add_index(labels, index_name=DNAStream.SNV, source_file=source_file)


    def batch_add_samples(self, labels, source_file=""):
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
        return self._batch_add_index(labels, index_name=DNAStream.SAMPLE, source_file=source_file)


    def _batch_add_index(self, labels, index_name, source_file=""):
        """
        Internal method to batch add multiple labels (SNVs or samples) to the dataset.

        Parameters
        ----------
        labels : list of str
            List of labels to add.
        index_name : str
            The index type, either DNAStream.SNV or DNAStream.SAMPLE.
        source_file : str, optional
            Path to the source file from which the labels are being added (default is an empty string).

        Returns
        -------
        dict
            A dictionary mapping labels to their respective indices.
        """
        index_dict = self._load_index(index_name)
        pre_size = len(index_dict)

        indices = []
        new = 0
        for lab in labels:
            if lab in index_dict:
                indices.append(index_dict[lab])
            else:
                next_idx = len(index_dict)
                indices.append(next_idx)
                index_dict[lab] = next_idx
                new += 1

        indices = np.array(indices)
        labels = np.array(labels)
        
        # Sort indices and labels concurrently for HDF5 fancy indexing
        labels, indices = labels[np.argsort(indices)].tolist(), indices[np.argsort(indices)].tolist()

        if self.file[f"{index_name}/label"].shape[0] != len(index_dict):
            if index_name == DNAStream.SNV:
                self._resize_all(m=len(index_dict))
            if index_name == DNAStream.SAMPLE:
                self._resize_all(n=len(index_dict))

        self.file[f"{index_name}/label"][indices] = labels
        self._save_index(index_dict, index_name)
        post_size = len(index_dict)

        if new > 0:
            self._update_index_log(index_name, new, pre_size, post_size, operation="add", source_file=source_file)

        if self.verbose:
            print(f"#{new} items added to {index_name} index")

        return index_dict

    def _update_index_log(self, index_name, num, pre_size, post_size, operation, source_file):
        """
        Log index updates (SNV or sample) in the HDF5 file.

        This method records modifications to the index, including the number of new entries, 
        operation type (e.g., "add"), and relevant metadata.

        Parameters
        ----------
        index_name : str
            The index type, either DNAStream.SNV or DNAStream.SAMPLE.
        num : int
            The number of new entries added to the index.
        pre_size : int
            The size of the index before the operation.
        post_size : int
            The size of the index after the operation.
        operation : str
            The operation performed (e.g., "add", "remove").
        source_file : str
            Path to the file from which the data was added.

        Notes
        -----
        - The log entry includes a timestamp, the operation type, the user, and the hostname.
        - The log dataset is resized dynamically to accommodate new entries.
        """
        source_file = str(pathlib.Path(source_file).resolve())
        log = self.file[f"{index_name}/log"]
        current_size = log.shape[0]
        new_size = current_size + 1
        log.resize((new_size,))
        timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S").encode("utf-8")
        user = getpass.getuser().encode("utf-8")
        hostname = socket.gethostname().encode("utf-8")
        log[current_size] = (timestamp_str, num, pre_size, post_size, operation, user, hostname, source_file)


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
        return self._get_data(dataset_name=f"{DNAStream.SNV}/data", indices=indices)


    def get_snv_log(self):
        """
        Retrieve the SNV index update log.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for SNV updates.
        """
        return self._get_data(dataset_name=f"{DNAStream.SNV}/log")


    def get_sample_log(self):
        """
        Retrieve the sample index update log.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing log entries for sample updates.
        """
        return self._get_data(dataset_name=f"{DNAStream.SAMPLE}/log")

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
        return self._get_data(dataset_name=f"{DNAStream.SAMPLE}/data", indices=indices)


    def _get_data(self, dataset_name, indices=None):
        """
        Internal method to retrieve structured data from the HDF5 file.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to retrieve (e.g., "SNV/data").
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
        columns = dataset.attrs['columns']

        if indices is not None:
            dataset = dataset[indices, :]
        else:
            dataset = dataset[:]

        df = pd.DataFrame(dataset, columns=columns)

        # Convert byte strings to UTF-8
        for col in df.select_dtypes(include=["object"]):
            if df[col].dtype == object:
                df[col] = df[col].apply(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x)

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
        self._add_data(indices, df, dataset_name=f"{DNAStream.SNV}/data", source_file=source_file)


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
        self._add_data(indices, df, dataset_name=f"{DNAStream.SAMPLE}/data", source_file=source_file)


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
            raise IndexError("Invalid indices passed to add_data method, check indices and try again")

        # Convert df to structured array
        structured_data = np.zeros(len(indices), dtype=dataset.dtype)
        for col in df.columns:
            structured_data[col] = df[col].to_numpy()

        # Update data at provided indices
        dataset[indices] = structured_data

        #log modifications 
        self._log_dataset_modification(dataset_name, operation="update", source_file=source_file)


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
            raise IndexError("Invalid indices passed to update_data method, check indices and try again")

        # Read affected rows from HDF5 into memory
        temp_data = dataset[indices]  # Read only required rows

        # Update only the specified column
        temp_data[col] = value

        # Write back only modified rows
        dataset[indices] = temp_data
        self._log_dataset_modification(dataset_name, operation="update", source_file=source_file)
   

    def add_maf_files(self, fnames, **kwargs)  :
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


    def add_maf_file(self, fname,
                                 missing_values= ["Unknown", "Na", "N/A", "na", "nan", 
                                                    "NaN", "NAN", "NONE", "None", "", "__UNKNOWN__"],
                                 required_cols =["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                                             "Reference_Allele", "Tumor_Seq_Allele2", "Entrez_Gene_Id"] ):
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

        #read MAF file and extract key info
        maf = pd.read_table(fname, low_memory=False) 
      
        missing_cols = set(required_cols) - set(maf.columns)

        maf.replace(missing_values, pd.NA, inplace=True)
        if missing_cols:
            maf = maf.reindex(columns=maf.columns.tolist() + list(missing_cols), fill_value=pd.NA)
        
        
        
        maf["label"] = maf.iloc[:, :4].astype(str).agg(":".join, axis=1)

        column_dict = {
            "label" : "label",
            "Chromosome" : "chrom",
            "Start_Position" : "pos",
            "End_Position" : "end_pos",
            "Reference_Allele" : "ref_allele",
            "Tumor_Seq_Allele2" : "alt_allele",
            "Hugo_Symbol" : "hugo",
            "Entrez_Gene_Id": "gene"
        }

        maf.rename(columns = column_dict, inplace=True)

    
        snv_labels = maf["label"].tolist() 


        try: 

            
                snv_idx = self.batch_add_snvs(snv_labels, source_file=fname)

                #sort dataframe according to the newly assigned indices in DNAStream
                maf.loc[:,"snv_idx"] = maf["label"].map(snv_idx) 
                maf = maf.sort_values("snv_idx")
                indices = maf["snv_idx"].tolist()
                snv_data = maf[[val for _, val in column_dict.items()]]
              
    
                self.add_snv_data(indices, snv_data, source_file=fname)


        except Exception as e:
       
            self.close()
            raise Exception(e)

    @staticmethod
    def _parse_file(fname, sep_word='tree', nskip=0, sep="\t"):
        """
        Parses a text file containing multiple edge lists of trees.
        @param fname: str filename to be parsed
        @param sep_word: str a word contained in the line separating distinct trees (default 'tree')
        @param nksip: int number of rows to skip before parsing
        """
        tree_list = []
        with open(fname, 'r+') as file:
            new_tree = None 
            for idx, line in enumerate(file):
                if idx < nskip:
                    continue
                if sep_word in line:
                    if new_tree:
                        tree_list.append(new_tree)
                    new_tree = []

                else:
                    edge = [ int(e) for e in line.strip().split(sep)]
                    new_tree.append((edge[0], edge[1])) 
            if new_tree:
                tree_list.append(new_tree)
   
        return tree_list, None

    def _add_tree(self,  edge_list,  tree_type,data=None, index=None):
       
        trees = self.file[f"{DNAStream.TREE}/{tree_type}_trees"]
        save = False
        if index is None:
            index = self._load_index(f"{DNAStream.TREE}/{tree_type}_trees/index_map")
            save = True
   
        numtrees = len(index)
        new_size = numtrees + 1

        for key in trees:
            if key != "index_map":
                trees[key].resize((new_size,))
        label = f"tree{numtrees}"
        index[f"tree{numtrees}"] =numtrees

        if save:
            self._save_index(index, f"{DNAStream.TREE}/{tree_type}_tree/index_map")
       
     
        trees["trees"][numtrees] =np.array(edge_list, dtype=EDGE_LIST_DTYPE)
        
        if data:
            data["label"] = label
            trees["data"][numtrees] = data  
    
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
                    print(f"#Warning! Entries already exist from {source_file} in '{dataset_name}'")
                  
                return False
            else:
                return True
       

    def add_trees_from_file(self, fname, tree_type="SNV", method="", safe=True):
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

            dataset_name = f"{DNAStream.TREE}/{tree_type}_trees/data"

            if safe and not self._check_safe(dataset_name, source_file):
                print("#Warning! Attempting overwrite in Safe mode, use safe=F, to force append trees.")
                return 

 

            tree_index = self._load_index(f"{DNAStream.TREE}/{tree_type}_trees/index_map")

            # Parse tree file according to method
            if method == "conipher":
                tree_list, scores = self._parse_file(fname)
            else:
                tree_list, scores = self._parse_file(fname, nskip=1, sep=" ")

            data_dtype = SCHEMA[DNAStream.TREE][f"{tree_type}_trees"]["data"]
          
            for i, edge_list in enumerate(tree_list):
                if scores is None:
                    dat = np.array([("", method, np.nan, i, source_file)], dtype=data_dtype)
                else:
                    dat = np.array([("", method, scores[i], i, source_file)], dtype=data_dtype)

                self._add_tree(edge_list, tree_type, data=dat, index=tree_index)

            self._save_index(tree_index, f"{DNAStream.TREE}/{tree_type}_trees/index_map")
            self._log_dataset_modification(f"{DNAStream.TREE}/{tree_type}_trees", operation="update", source_file=fname)

            if self.verbose:
                print(f"#{len(tree_list)} {tree_type} trees added from {method}.")

        except Exception as e:
            self.close()
            raise Exception(e)


    def add_snv_tree_from_edge_list(self, edge_list, method="", score=np.nan, rank=np.nan):
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
        data_dtype = SCHEMA[DNAStream.TREE]["SNV_trees"]["data"]
        dat = np.array([("", method, score, rank, "")], dtype=data_dtype)
        self._add_tree(edge_list, "SNV", data=dat)
        self._log_dataset_modification(f"{DNAStream.TREE}/SNV_tree", operation="update")

    def add_pyclone_file(self, fname):
        pass 



    def close(self):
        """
        Close the HDF5 file safely.

        This should be called when finished working with the DNAStream object.
        """
      
        self.file.close()
        if self.verbose:
            print(f"#Stream to connection {self.filename} closed.")


