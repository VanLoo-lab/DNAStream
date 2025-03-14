import json
import pathlib
import getpass
import socket
from datetime import datetime
import h5py
import numpy as np


class BaseIndex:
    """Base class for handling index storage in an HDF5 dataset."""

    def __init__(self, file, name, metadata_dtype="", verbose=False):
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

        # Ensure the index dataset exists
        if self.label_name not in self.file:
            self.file.create_dataset(
                self.label_name,
                shape=(0,),
                maxshape=(None,),
                dtype=h5py.string_dtype("utf-8"),
            )

        if self.dat_name not in self.file:
            self.file.create_dataset(
                self.dat_name,
                shape=(0,),
                maxshape=(None,),
                dtype=self.dat_dtype,
            )
            self.file[self.dat_name].attrs["columns"] = self.dat_dtype.names

        self.labels = self.file[self.label_name]
        self.metadata = self.file[self.dat_name]

        # Load index into memory for fast access
        self._labels_cache = list(self.labels[:])  # Convert NumPy array to list
        self._index_cache = {lab: i for i, lab in enumerate(self._labels_cache)}
        self._metadata_cache = self.metadata[:]  # numpy array

    def size(self):
        """Return the number of labels in the index."""
        return len(self._labels_cache)

    def _resize(self, new_size):
        """Resize the index dataset."""
        self.labels.resize((new_size,))
        self.metadata.resize((new_size,))

    def _save_index(self):
        """Save index to HDF5 dataset if modified."""
        if self._modified:
            self._resize(new_size=self.size())

            # write the cache to disk
            self.labels[:] = self._labels_cache  # Update HDF5 dataset
            self.metadata[:] = self._metadata_cache  # Expand dataset
            self._modified = False  # Reset modification flag

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

    def add(self, labels, metadata={}):
        """
        Add multiple labels to the index.

        Parameters
        ----------
        labels : list of str
            List of new labels to be added.
        """
        new_idx = self.size()
        added_indices = {}
        for label in labels:
            if label not in self._index_cache:  # Avoid duplicates
                self._labels_cache.append(label)
                self._index_cache[label] = new_idx
                added_indices[label] = new_idx
                new_idx += 1
                self._modified = True  # Mark dataset as modified

        new_metadata = np.zeros(len(added_indices), dtype=self.dat_dtype)

        # Assign metadata values (use default if not provided)
        for i, (label, idx) in enumerate(added_indices.items()):
            entry_metadata = metadata.get(
                label, {}
            )  # Get metadata dict for label or empty dict
            for field in self.dat_dtype.names:  # Loop through all fields in dtype
                new_metadata[i][field] = entry_metadata.get(
                    field, self._default_metadata_value(field)
                )

        # Insert new metadata into the cache
        self._metadata_cache = np.concatenate((self._metadata_cache, new_metadata))

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

    def get_metadata(self, label):
        """Get metadata for a given label.

        Parameters
        ----------
        label : str
            The label to look up.

        """
        idx = self._index_cache.get(label, None)
        return self.metadata[idx] if idx is not None else None


class GlobalIndex(BaseIndex):
    """Handles global indices (SNV, sample)."""

    def __init__(self, file, name, metadata_dtype, verbose=False):
        super().__init__(file, name, metadata_dtype, verbose)
        self.log_name = f"{self.group}/log"
        self.cluster_name = f"{self.group}/cluster"

        # Check if datasets exist, create if not
        if self.log_name not in self.file:
            self.file.create_dataset(
                self.log_name, (0,), maxshape=(None,), dtype=h5py.string_dtype("utf-8")
            )
            self.file[self.log_name].attrs["columns"] = self.dat_dtype.names

        if self.cluster_name not in self.file:
            self.file.create_dataset(
                self.cluster_name, (0,), maxshape=(None,), dtype="i"
            )

        self.log = self.file[self.log_name]
        self.cluster = self.file[self.cluster_name]

    def get_log(self):
        """Return the log dataset."""
        return self.log

    def _resize(self, new_size):
        super()._resize(new_size)
        self.cluster.resize((new_size,))

    def update_cluster(self, label, cluster):
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

    def update_clusters(self, cluster_dict):
        """
        Update the cluster assignments for a list of labels.

        Parameters
        ----------
        cluster_dict
        """
        for label, cluster in cluster_dict.items():
            self.update_cluster(label, cluster)

    def add(self, labels, metadata=None, clusters=None, source_file=""):
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
        super().add(labels, metadata)
        post_size = self.size()

        if clusters:
            if len(clusters) != len(labels):
                raise ValueError("len(clusters) must match number of len(labels).")

            cluster_dict = dict(zip(labels, clusters))
            self.update_cluster(cluster_dict)

        if post_size > pre_size:
            self._update_index_log(
                post_size - pre_size, pre_size, post_size, "add", source_file
            )

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

        log_entry = f"{timestamp_str},{num},{pre_size},{post_size},{operation},{user},{hostname},{source_file}"
        self.log[current_size] = log_entry

    def save_index(self):
        """Force save the index to disk."""
        self._modified = True
        self._save_index()
