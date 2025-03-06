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

#TODO
# - create mkdocs
# - add pyclone addition
# - add bulk phylogenies as edge lists

"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing quality scores, number of callers, etc
     |── cluster
     |── index_map             #json string for fast loading and saving
     |-- log
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
     |-- log
 ├── read_counts/               # Read count matrices
 │   ├── bulk/                  # Bulk sequencing read counts
 │   │   ├── variant       # SNVs x Samples (variant read counts)
 │   │   ├── total         # SNVs x Samples (total read counts)
 │   ├── lcm/                    # LCM sequencing read counts
 │   │   ├── variant       
 │   │   ├── total         
 │   ├── scdna/                  # scDNA-seq read counts
 │   │   ├── variant       
 │   │   ├── total         
 ├── metadata/                   # Metadata storage
 │   ├── sample_info              # Sample IDs
 │   ├── processing_parameters
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

    Class Attributes
    ----------------
    BULK : str
        Label for bulk sequencing data ("bulk").
    LCM : str
        Label for laser capture microdissection (LCM) sequencing data ("lcm").
    SCDNA : str
        Label for single-cell DNA sequencing data ("scdna").
    SNV : str
        Label for the SNV (single-nucleotide variant) index ("SNV").
    SAMPLE : str
        Label for the sample index ("sample").
    META_TABLES : list
        List of metadata tables stored within each index. Includes:
        - "data" : Contains structured SNV/sample metadata.
        - "label" : Stores unique labels for SNVs and samples.
        - "cluster" : Stores cluster assignments for SNVs or samples.
        - "index_map" : Maps SNV/sample labels to indices.
        - "log" : Stores operation logs for tracking modifications.
    MODALITIES : list
        List of supported sequencing modalities: ["bulk", "lcm", "scdna"].
    INDICES : list
        List of available indices: ["SNV", "sample"].

    Instance Attributes
    -------------------
    filename : str
        Path to the HDF5 file used for data storage.
    verbose : bool
        If True, enables verbose output during operations (default: False).
    file : h5py.File
        The HDF5 file object that manages storage and retrieval.
    schema : dict
        A dictionary defining the structure and data types of SNV and sample datasets.

    Methods
    -------
    __init__(filename, verbose=False, snv_dtype, sample_dtype)
        Initializes the HDF5 storage and creates necessary datasets if they do not exist.
    add_read_counts(fname, source, location=None)
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

    META_TABLES = ["data", "label", "cluster", "index_map", "log"]
    MODALITIES = [BULK, LCM, SCDNA]
    INDICES = [SNV, SAMPLE]

    def __init__(self, filename, verbose=False, 
            snv_dtype = np.dtype([
                    ("label", h5py.string_dtype(encoding="utf-8")),
                    ("chrom", "S5"),  # Fixed-length string (5 characters max)
                    ("pos", "i8"), 
                    ("end_pos", "i8"),    # Integer
                    ("ref_allele", "S1"),  # Single-character string
                    ("alt_allele", "S1"),  # Single-character string
                    ("hugo", "S15"), # Fixed-length string (15 characters max)
                    ("gene", "S10"), #Fixed-length string (10 characters ma)       
             ]),
            sample_dtype = np.dtype([
                    ("label", h5py.string_dtype(encoding="utf-8")),  #sample description
                    ("patient", "S10"),  # Fixed-length string (10 characters max)
                    ("source", "S10"),  # Fixed-length string (10 characters max)
                    ("location", "S15") , # Fixed-length string (15 characters max)
                    ("bam_file", h5py.string_dtype(encoding="utf-8")),  # Single-character string

            ])
    ):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose

    

        log_dtype = np.dtype([
            ("timestamp", h5py.string_dtype(encoding="utf-8")),
            ("number", "i8"),
            ("index_size_before", "i8"),
            ("index_size_after", "i8"),
            ("operation", "S15"),
            ("user", "S15"),
            ("hostname", "S15"),
            ("file",  h5py.string_dtype(encoding="utf-8"))
        ])

        


        self.file = h5py.File(filename, "a")  # Append mode (does not overwrite)
        

        # snv_grp = self.self.h5file.create_group(SNV)
        # snv_grp = snv_grp.attrs['label'] = np.array(shape=0, dtype=str)
        # snv_grp = snv_grp.attrs['cluster']= np.array(shape=0, dtype=np.int32)

        self.schema = {DNAStream.SNV: 
                {"label" :  h5py.string_dtype(encoding="utf-8"),
                "cluster" : "i8",
                "data" : snv_dtype,
                "index_map": h5py.string_dtype("utf-8"),
                "log" : log_dtype,
                },
            DNAStream.SAMPLE : 
                {
                "label":    h5py.string_dtype(encoding="utf-8") ,
                "cluster" : "i8",
                "data" : sample_dtype,
                "index_map":    h5py.string_dtype("utf-8"),
                "log" : log_dtype
                }
        }

        # Create shared SNV index if not already present
        for index in DNAStream.INDICES:
            for data in DNAStream.META_TABLES:
                if f"{index}/{data}" not in self.file:
                    self.file.create_dataset(f"{index}/{data}", shape=(0,), maxshape=(None,), 
                                        dtype=self.schema[index][data], compression="gzip")
                    if data in ["data", "log"]:   
                        self.file[f"{index}/{data}"].attrs['columns'] = list(self.schema[index][data].names)

        # Create group structure
        for modality in DNAStream.MODALITIES:
                group_path = f"read_counts/{modality}"
                if group_path not in self.file:
                    self.file.create_group(group_path)
                for reads in ["variant", "total"]:
                    if f"{group_path}/{reads}" not in self.file:
                        self.file.create_dataset(f"{group_path}/{reads}", shape=(0,0), maxshape=(None,None), dtype='i', 
                                                 compression="gzip", chunks=(1, 5000))
                        self.file[f"{group_path}/{reads}"].dims[0].label = DNAStream.SNV
                        self.file[f"{group_path}/{reads}"].dims[1].label = DNAStream.SAMPLE


        
    def __str__(self):
        """To string method"""
        m = self.file[f"{DNAStream.SNV}/label"].shape[0]
        n = self.file[f"{DNAStream.SAMPLE}/label"].shape[0]

        mystr = f"DNAStream object with {m} SNVs and {n} samples" 
        mystr += f"\nHDF5 File: {self.filename}"
      
        return mystr
    


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
            for snv_data in DNAStream.META_TABLES:
                if snv_data == "log":
                    continue
                group_path = f"{DNAStream.SNV}/{snv_data}"
            
                self.file[group_path].resize((m,))
        else:
            m= self.file[f"{DNAStream.SNV}/label"].shape[0]

        if n:
            for sample_data in DNAStream.META_TABLES:
                if sample_data == "log":
                    continue
                group_path = f"{DNAStream.SAMPLE}/{sample_data}"
                self.file[group_path].resize((n,))
        else:
            n = self.file[f"{DNAStream.SAMPLE}/label"].shape[0]


        for modality in DNAStream.MODALITIES:
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
        return self._save_index(index_dict, DNAStream.SNV)


    def save_sample_index(self, index_dict):
        """
        Save the sample index to the HDF5 file.

        Parameters
        ----------
        index_dict : dict
            Dictionary mapping sample labels to their respective indices.
        """
        return self._save_index(index_dict, DNAStream.SAMPLE)


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
        if f"{index_name}/index_map" in self.file:
            index_data = self.file[f"{index_name}/index_map"][()]
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
        self.file[f"{index_name}/index_map"].resize((1,))  # Ensure space in dataset
        self.file[f"{index_name}/index_map"][0] = index_json  # Store JSON string in HDF5


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


    def add_snv_data(self, indices, df):
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
        self._add_data(indices, df, dataset_name=f"{DNAStream.SNV}/data")


    def add_sample_data(self, indices, df):
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
        self._add_data(indices, df, dataset_name=f"{DNAStream.SAMPLE}/data")


    def _add_data(self, indices, df, dataset_name):
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

        # Convert DataFrame to structured NumPy array
        structured_data = np.zeros(len(indices), dtype=dataset.dtype)
        for col in df.columns:
            structured_data[col] = df[col].to_numpy()

        # Update data at provided indices
        dataset[indices] = structured_data


    def _update_value(self, indices, dataset_name, col, value):
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
              
    
                self.add_snv_data(indices, snv_data)


        except Exception as e:
       
            self.close()
            raise Exception(e)



    
    def close(self):
        """
        Close the HDF5 file safely.

        This should be called when finished working with the DNAStream object.
        """
        self.file.close()


############ TEST ###################

# rc_pth = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/scdna/read_counts"
# ds = DNAStream(filename="temp.h5", verbose=True)
# print(ds)
# samples = [f"GEM2.2_PT_{i}" for i in range(2,6)] 
# maf_files = [f"../data_summary/WGS_MUTATION/Somatic/2outof3_SNV/{sample}_SNVs_2outof3.maf" for sample in samples]
# indices = ds.add_maf_files(maf_files)
# print(ds.get_snv_log())
# print(ds)
# read_count_file = f"{rc_pth}/GEM2.2.csv"
# ds.add_read_counts(read_count_file, source="scdna")
# print(ds.get_snv_log())
# print(ds.get_sample_log())
# print(ds)
# ds.close()
