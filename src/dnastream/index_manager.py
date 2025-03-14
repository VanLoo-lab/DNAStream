import json
import pathlib
import getpass
import socket
from datetime import datetime
import h5py

class BaseIndex:
    """Base class for handling index storage in an HDF5 dataset."""

    def __init__(self, file, name, verbose=False):
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
        self.file = file
        self.table = name
        self.verbose = verbose
        self.name = f"{self.table}/labels"
        self._modified = False

        # Ensure the index dataset exists
        if self.name not in self.file:
            self.file.create_dataset(
                self.name, shape=(0,), maxshape=(None,), dtype=h5py.string_dtype("utf-8")
            )

        self.labels = self.file[self.name]
        
        # Load index into memory for fast access
        self._labels_cache = list(self.labels[:])  # Convert NumPy array to list
        self._index_cache = {lab: i for i, lab in enumerate(self._labels_cache)}

    def size(self):
        """Return the number of labels in the index."""
        return len(self._labels_cache)
    
    def _save_index(self):
        """Save index to HDF5 dataset if modified."""
        if self._modified:
            new_size = self.size()
            self.file[self.name].resize((new_size,))  # Expand dataset
            self.file[self.name][:] = self._labels_cache  # Update HDF5 dataset
            self._modified = False  # Reset modification flag

    def add_labels(self, labels):
        """
        Add multiple labels to the index.

        Parameters
        ----------
        labels : list of str
            List of new labels to be added.
        """
        new_idx = self.size()
        for label in labels:
            if label not in self._index_cache:  # Avoid duplicates
                self._labels_cache.append(label)
                self._index_cache[label] = new_idx
                new_idx += 1
                self._modified = True  # Mark dataset as modified

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

    def _allocate_labels(self, num, prefix=None ):
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
            prefix = self.table
        new_idx = self.size()
        labels = [f"{prefix}_{i}" for i in range(new_idx, new_idx + num)]
        return self.add_labels(labels)
    
class GlobalIndex(BaseIndex):
    """Handles global indices (SNV, sample)."""
    
    def __init__(self, file, name, verbose=False):
        super().__init__(file, name, verbose)
        self.log_name = f"{self.table}/log"

        # Check if datasets exist, create if not
        if self.log_name not in self.file:
            self.file.create_dataset(self.log_name, (0,), maxshape=(None,), dtype=h5py.string_dtype("utf-8"))

        self.log = self.file[self.log_name]

    def get_log(self):
        """Return the log dataset."""
        return self.log

    def add_labels(self, labels, source_file=""):
        """
        Add multiple labels and log the operation.

        Parameters
        ----------
        labels : list of str
            List of labels to be added.
        source_file : str, optional
            The source file that triggered the modification.
        """
        pre_size = self.size()
        super().add_labels(labels)
        post_size = self.size()

        if post_size > pre_size:
            self._update_index_log(post_size- pre_size, pre_size, post_size, "add", source_file)

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