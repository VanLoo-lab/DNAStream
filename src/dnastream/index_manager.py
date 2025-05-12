import json
import time
import pathlib
import getpass
import socket
from datetime import datetime
import h5py
import numpy as np
import pandas as pd
from .datatypes import LOG_DTYPE


class BaseIndex:
    """Base class for handling index storage in an HDF5 dataset."""

    def __init__(self, file, name, metadata_dtype, tracked_tables=None, verbose=False):
        """
        Initialize an index object.

        Parameters
        ----------
        file : h5py.File
            The HDF5 file where the index dataset is stored.
        name : str
            The index name (e.g., 'SNV', 'sample', 'trees', 'copy_number').
        verbose : bool, optional
            Whether to print debug messages (default is False).
        """

        # if "index" not in file:
        #     raise ValueError(f"index group not found in file!")

        self.file = file
        self.group = name
        self.verbose = verbose
        self.label_name = f"{self.group}/labels"
        self.dat_name = f"{self.group}/metadata"
        self.dat_dtype = metadata_dtype
        self._modified = False
        self.last_saved_timestamp = None

        # Ensure the index dataset exists
        if self.label_name not in self.file:
            self.file.create_dataset(
                self.label_name,
                shape=(0,),
                maxshape=(None,),
                dtype=h5py.string_dtype("utf-8"),
                chunks=True,
            )

        if self.dat_name not in self.file:

            self.file.create_dataset(
                self.dat_name,
                shape=(0,),
                maxshape=(None,),
                dtype=self.dat_dtype,
                chunks=True,
            )
            self.file[self.dat_name].attrs["columns"] = self.dat_dtype.names

        self.labels = self.file[self.label_name]
        self.metadata = self.file[self.dat_name]

        # Load index into memory for fast access
        self._labels_cache = [lab.decode("utf-8") if isinstance(lab, bytes) else lab for lab in self.labels[:]] # Convert NumPy array to list
        # print(f"SIZE OF LABELS CACHE: {len(self._labels_cache)}")
        self._index_cache = {lab: i for i, lab in enumerate(self._labels_cache)}
        self._metadata_cache = self.metadata[:]  # numpy array

        # Initialize tracked tables
        self.tracked_tables = {}
        if tracked_tables:
            for table_name, axis in tracked_tables:
                self.add_tracked_table(table_name, axis)

    def size(self):
        """Return the number of labels in the index."""
        return len(self._labels_cache)

    def label_to_idx(self, label):
        """Return the index of a label."""
        return self._index_cache.get(label, None)

    def _resize(self, new_size):
        """Resize the index dataset."""
        self.labels.resize((new_size,))
        self.metadata.resize((new_size,))

    def flush_metadata(self):
        """Write metadata cache to disk."""
        if self.metadata.shape[0] != len(self._metadata_cache):
            self.metadata.resize((len(self._metadata_cache),))
        self.metadata[:] = self._metadata_cache

    def flush_labels(self):
        """Write metadata cache to disk."""
        self._resize(new_size=self.size())
        # write the cache to disk
        self.labels[:] = self._labels_cache

    def _save_index(self):
        """Save index to HDF5 dataset if modified."""
        if self._modified:
            self.flush_labels()

            # Update HDF5 dataset

            self.flush_metadata()
            # self.metadata[:] = self._metadata_cache  # Expand dataset
            self._modified = False  # Reset modification flag

            # Update last saved timestamp
            self.last_saved_timestamp = time.time()

            if self.verbose:
                print(
                    f"Index saved at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.last_saved_timestamp))}"
                )

    def _default_metadata_value(self, field):
        """Returns a default value for missing metadata fields based on dtype."""
        field_type = self.metadata.dtype[field]

        if np.issubdtype(field_type, np.integer):
            return -1  # Default integer value
        elif np.issubdtype(field_type, np.floating):
            return np.nan  # Default float value
        elif np.issubdtype(field_type, np.str_) or field_type == object:
            return ""  # Use an empty string for variable-length UTF-8 strings
        elif np.issubdtype(field_type, np.bytes_):
            return b""  # Use bytes for fixed-length strings
        else:
            raise ValueError(f"Unhandled metadata field type: {field_type}")

    def _check_and_convert_dtype(self, field, value):
        """Check and convert value to match expected dtype for metadata insertion."""
        dtype = self.dat_dtype[field]
        try:
            if isinstance(dtype, h5py.Datatype) or dtype.kind == "O":
                return str(value)
            elif np.issubdtype(dtype, np.bytes_):
                return str(value).encode("utf-8")
            elif np.issubdtype(dtype, np.integer):
                return int(value)
            elif np.issubdtype(dtype, np.floating):
                return float(value)
            else:
                return value
        except Exception as e:
            raise ValueError(
                f"Failed to convert field '{field}' with value '{value}' to dtype '{dtype}': {e}"
            )

    def add_tracked_table(self, table_name, axis):
        if table_name not in self.file:
            raise ValueError(f"Dataset {table_name} must already exists in file.")

        if self.file[table_name].shape[axis] != self.size():
            raise ValueError(
                f"Dataset {table_name} shape for axis {axis} does not match the size of the current index."
            )

        self.tracked_tables[table_name] = axis

    def _get_data(self, dataset, indices=None):
        """
        Internal method to retrieve structured data from the HDF5 file.

        Parameters
        ----------
        dataset : reference to a HDF5 dataset
            The name of the dataset to retrieve (e.g., "SNV/metadata").

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the requested dataset.

        Notes
        -----
        - The method decodes byte strings into UTF-8.
        """

        columns = dataset.attrs["columns"]
        if indices:
            dataset = dataset[indices]
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

    def insert_metadata(self, metadata_dict):
        """
        Insert or update metadata for existing labels.

        Parameters
        ----------
        metadata_dict : Dict[str, Dict[str, Any]]
            Mapping from labels to a dict of metadata field names and values.
            Missing fields are filled with default values during index creation.
            Fields not included in a metadata row will be left unchanged.
        """
        if not metadata_dict:
            return

        # Ensure all labels are in the index (adds them if not already there)
        added_indices = self.add(list(metadata_dict.keys()))  # {label: index}

        # For each label, update metadata in-place
        for label, idx in added_indices.items():
            row_data = metadata_dict[label]
            for field, value in row_data.items():
                if field in self.dat_dtype.names:

                    self._metadata_cache[idx][field] = self._check_and_convert_dtype(
                        field, value
                    )
                else:
                    if self.verbose:
                        print(
                            f"#Warning: field '{field}' not in metadata dtype for {self.group}. Skipping."
                        )
        self._modified = True

    def add(self, labels):
        """
        Add multiple labels to the index.

        Parameters
        ----------
        labels : list of str
            List of new labels to be added.
        """
        pre_size = self.size()
        post_size = self.size()
        added_indices = {}
        for label in labels:
            if label not in self._index_cache:  # Avoid duplicates
                self._labels_cache.append(label)
                self._index_cache[label] = post_size
                added_indices[label] = post_size
                post_size += 1
                self._modified = True  # Mark dataset as modified

        if pre_size != post_size:
            self._resize_tracked_tables(post_size)

        for _ in range(post_size - pre_size):
            default_row = tuple(
                self._default_metadata_value(f) for f in self.dat_dtype.names
            )
            self._metadata_cache = np.append(
                self._metadata_cache, np.array([default_row], dtype=self.dat_dtype)
            )

        self._save_index()  # Save if modified
        return {label: self._index_cache[label] for label in labels}

    def get_labels(self):
        """Return the list of labels."""
        return self._labels_cache

    def get_index(self):
        """Return the current index dictionary."""
        return self._index_cache  # Use the cached version

    ### **Dictionary-Like Behavior**
    def __getitem__(self, key):
        """Enable dictionary-style lookups."""
        return self._index_cache.get(key, None)

    def __contains__(self, key):
        """Enable 'in' keyword for checking membership."""
        return key in self._index_cache

    def to_json(self):
        """Return the index as a JSON string."""
        return json.dumps(self._index_cache, indent=4)

    def indices_by_label(self, labels):
        """
        Retrieve indices for a list of labels.

        Parameters
        ----------
        labels : list of str
            List of labels to look up.

        Returns
        -------
        list of int
            List of indices corresponding to the labels. None if labels are not in index.
        """
        return [self._index_cache.get(label, None) for label in labels]

    def resize(self, num, prefix=None):
        """
        Allocate a block of labels.

        Parameters
        ----------
        num : int
            The number of labels to allocate.

        Returns
        -------
        list of int
            List of indices for the allocated labels.
        """
        if not prefix:
            prefix = self.group
        new_idx = self.size()
        labels = [f"{prefix}_{i}" for i in range(new_idx, new_idx + num)]
        return self.add(labels)

    def _resize_tracked_tables(self, new_size):
        """Resize all tables tracked by this index when its size changes."""

        for table_name, axis in self.tracked_tables.items():
            if table_name in self.file:
                shape = list(self.file[table_name].shape)
                if shape[axis] != new_size:
                    shape[axis] = new_size  # Resize only the correct axis
                    self.file[table_name].resize(tuple(shape))
            else:
                if self.verbose:
                    print(
                        f"Warning {table_name} not in file! skipping resizing. Modify tracked tables for {self.group} index to suppress this warning."
                    )

    def get_metadata(self, labels=None, format="numpy"):
        """
        Get metadata for specified labels in a given format.

        Parameters
        ----------
        labels : list of str, optional
            Labels to retrieve metadata for. If None, returns all.
        format : str, optional
            Output format: "numpy", "dict", or "pandas". Default is "numpy".

        Returns
        -------
        np.ndarray, dict, or pd.DataFrame
            Metadata in the requested format.
        """
        if labels:
            idxs = [self[label] for label in labels]
            data = self._metadata_cache[idxs]
            label_subset = labels
        else:
            data = self._metadata_cache
            label_subset = self.get_labels()

        if format == "numpy":
            return data
        elif format == "dict":
            return {
                label: {field: row[field] for field in data.dtype.names}
                for label, row in zip(label_subset, data)
            }
        elif format == "pandas":
            df = pd.DataFrame(data)
            for col in df.select_dtypes(include=["object"]):
                df[col] = df[col].apply(
                    lambda x: x.decode("utf-8") if isinstance(x, bytes) else x
                )
            return df
        else:
            raise ValueError(f"Unsupported format: {format}")


class GlobalIndex(BaseIndex):
    """Handles global indices (SNV, sample)."""

    def __init__(self, file, name, metadata_dtype, tracked_tables, verbose=False):
        super().__init__(file, name, metadata_dtype, tracked_tables, verbose)
        self.log_name = f"{self.group}/log"
        self.cluster_name = f"{self.group}/cluster"
        self.log_header = [
            "Timestamp",
            "Entries Added",
            "Pre-size",
            "Post-size",
            "Operation",
            "User",
            "Host",
            "Source File",
        ]

        # Check if datasets exist, create if not
        if self.log_name not in self.file:
            self.file.create_dataset(
                self.log_name, (0,), maxshape=(None,), dtype=LOG_DTYPE
            )
            self.file[self.log_name].attrs["columns"] = self.log_header

        if self.cluster_name not in self.file:
            self.file.create_dataset(
                self.cluster_name, (0,), maxshape=(None,), dtype="i"
            )

        self.log = self.file[self.log_name]
        self._update_index_log(0, 0, 0, "create")
        self.cluster = self.file[self.cluster_name]

    def get_log(self):
        """Return the log dataset."""
        return self._get_data(self.log)

    def _resize(self, new_size):
        super()._resize(new_size)
        self.cluster.resize((new_size,))

    def _update_cluster(self, label, cluster):
        """
        Update the cluster assignment for a label.

        Parameters
        ----------
        label : str
            The label to update.
        cluster : int
            The cluster assignment.
        """
        if label in self:
            self.cluster[self[label]] = cluster
        else:
            if self.verbose:
                print(f"Warning! Label {label} not found in index, skipping cluster.")

    def update_clusters(self, cluster_dict, source_file=""):
        """
        Update the cluster assignments for a list of labels.

        Parameters
        ----------
        cluster_dict
        """
        for label, cluster in cluster_dict.items():
            self._update_cluster(label, cluster)

        self._update_index_log(
            len(cluster_dict),
            self.size(),
            self.size(),
            "modified clusters",
            source_file,
        )

    def get_clusters(self):
        """Return the cluster assignments."""

        return {label: self.cluster[self[label]] for label in self._labels_cache}

    def add(self, labels, source_file=""):
        """
        Add multiple labels and log the operation.

        Parameters
        ----------
        labels : list of str
            List of labels to be added.
        metadata : dict, optional
            Metadata to be added for each label.
        clusters : list of int, optional
            Cluster assignments for each label.
        source_file : str, optional
            The source file that triggered the modification.
        """
        pre_size = self.size()
        index_dict = super().add(labels)
        post_size = self.size()

        # if clusters:
        #     if len(clusters) != len(labels):
        #         raise ValueError("len(clusters) must match number of len(labels).")

        #     cluster_dict = dict(zip(labels, clusters))
        #     self._update_cluster(cluster_dict)

        if post_size > pre_size:
            # Resize dependent datasets
            self._update_index_log(
                post_size - pre_size, pre_size, post_size, "add", source_file
            )
        return index_dict

    def _update_index_log(self, num, pre_size, post_size, operation, source_file=""):
        """
        Log index updates (SNV or sample) in the HDF5 file.

        Parameters
        ----------
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

        current_size = self.log.shape[0]
        new_size = current_size + 1
        self.log.resize((new_size,))
        timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        user = getpass.getuser()
        hostname = socket.gethostname()

        log_entry = np.array(
            [
                (
                    timestamp_str,
                    num,
                    pre_size,
                    post_size,
                    operation,
                    user,
                    hostname,
                    source_file,
                )
            ],
            dtype=LOG_DTYPE,
        )
        self.log[current_size] = log_entry

    def save_index(self):
        """Force save the index to disk."""
        self._modified = True
        self._save_index()

    def print_log(self):
        """
        Print the index log in a readable table format.
        """
        if self.log.shape[0] == 0:
            print("Log is empty.")
            return

        # Extract field names dynamically from dtype
        headers = self.log.dtype.names
        col_widths = [
            max(len(h), 12) for h in headers
        ]  # Ensure minimum width for readability

        # Print headers
        header_str = " | ".join(h.ljust(w) for h, w in zip(headers, col_widths))
        print(header_str)
        print("-" * len(header_str))  # Separator line

        # Print log entries
        for entry in self.log:
            values = [
                str(entry[h]) for h in headers
            ]  # Convert structured array row to list of strings
            formatted_entry = " | ".join(v.ljust(w) for v, w in zip(values, col_widths))
            print(formatted_entry)


class LocalIndex(BaseIndex):
    """Handles local indices (tree structures, copy numbers, etc.)."""

    def __init__(self, file, name, metadata_dtype, tracked_tables, verbose=False):
        super().__init__(file, name, metadata_dtype, tracked_tables, verbose)

    def get_metadata(self):
        """Return metadata for this local index."""
        return self._get_data(self.metadata)

    def get_labels(self):
        """Return labels for this local index."""
        return self._labels_cache


class DependentIndexView:
    def __init__(self, global_index, tracked_tables,  predicate_fn, file, idx_view_mapping_table,verbose=False):
        self.global_index = global_index
        self.tracked_tables = tracked_tables  # list of (table_name, axis)
        self.predicate_fn = predicate_fn
        self.file = file
        self.verbose = verbose
        self.idx_view_mapping_table = idx_view_mapping_table
        if self.idx_view_mapping_table in self.file:
                self._local_to_global_idx = self.file[self.idx_view_mapping_table][:].tolist()
         
        else:
            raise ValueError(f"Error! {idx_view_mapping_table} is not in datafile.")
   

        self.labels = [self.global_index._labels_cache[idx] for idx in self._local_to_global_idx]

  
    def label_to_view_idx(self,label):
        return self.labels.index(label)

    def label_to_global_idx(self, label):
        return self.global_index.label_to_idx(label)
    

    def register(self, labels):
        """Add labels that match the predicate and resize only relevant tracked tables.
        DOES NOT ADD LABELS TO GLOBAL IDX!!!
        """
        metadata = self.global_index.get_metadata(labels)
        filtered_labels = [
            l for l, row in zip(labels, metadata) if self.predicate_fn(row)
        ]

  
        modified = False
        for f in filtered_labels:
            if f not in self.labels:
                modified =True
                self.labels.append(f)
                self._local_to_global_idx.append(self.label_to_global_idx(f))
        
        # if self.verbose:
        #     print(
        #         f"#Filtering {len(filtered_labels)} matching labels for DependentIndexView..."
        #     )
        if modified:
            ds = self.file[self.idx_view_mapping_table]
            if ds.shape[0] != len(self._local_to_global_idx):
                ds.resize((len(self._local_to_global_idx),))
            ds[:] = np.array(self._local_to_global_idx)

            for table_name, axis in self.tracked_tables:
                if table_name not in self.file:
                    continue
                current_shape = self.file[table_name].shape
                if axis == 0:
                    new_shape = (current_shape[0] + len(filtered_labels), current_shape[1])
                else:
                    new_shape = (current_shape[0], current_shape[1] + len(filtered_labels))
                if current_shape != new_shape:
                    self.file[table_name].resize(new_shape)

        return {f : self.label_to_global_idx(f) for f in filtered_labels}
