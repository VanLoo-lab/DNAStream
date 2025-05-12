
import os
import glob
import pathlib
import re
import csv
import pandas as pd
import h5py
import numpy as np
from enum import Enum


# import json
from .index_manager import LocalIndex, GlobalIndex, DependentIndexView
from .enums import GlobalIndexName, LocalIndexName, Modalities, TreeType, SchemaGroups


from .datatypes import EDGE_LIST_DTYPE, ALLELE_SPECIFIC_CN_DTYPE

from .utils import (
    wrap_list,
    require_file_exists

)

from .schema import (
    SCHEMA,
    STRUCT_ARRAYS,
    LOCAL_INDEX,
    GLOBAL_INDEX,
    COPY_NUMBER_LAYER_DICT,
    get_schema_value,
)

class IOMixin:
    """
    Mixin class that provides I/O operations for reading and writing various types
    of genomic data into the DNAStream HDF5 structure.

    This class includes methods for:
    - Loading SNV, copy number, and tree data from common file formats (e.g., MAF, ASCAT, Battenberg, PyClone)
    - Inserting read counts and sample metadata
    - Parsing input directories and matching filenames via regex
    - Safely handling errors and file inconsistencies

    The IOMixin is intended to be used as part of the DNAStream class, not standalone.
    """

    # -----------------------------------
    # Index and metadata I/O functions
    # -----------------------------------
    def add_pyclone_clusters(self, source_file):
        """
        Parse a PyClone standard output file and extract the SNV cluster assignments.
        Updates the SNV global index cluster assignments accordingly.

        Parameters
        ----------
        source_file : str
            Path to the PyClone file to parse.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If the PyClone file cannot be parsed or required columns are missing.

        Notes
        -----
        The PyClone results file is expected to be tab-delimited and include columns:
        mutation_id, sample_id, cluster_id, cellular_prevalence, cellular_prevalence_std, cluster_assignment_prob.
        """

        # TODO: What do do about CCFs and cluster assignment probs?

        # if not self.in_context:
        #     raise RuntimeError("DNAStream I/O functions must be used inside a context manager.")
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

    def add_maf_files(self, fnames, **kwargs):
        """
        Add multiple MAF (Mutation Annotation Format) files to the SNV index and update the index log.

        This method iterates over a list of MAF file paths and processes each file
        using `add_maf_file`, updating the SNV index and data in the HDF5 dataset.

        Parameters
        ----------
        fnames : list of str
            A list of file paths to MAF files.
        **kwargs : dict, optional
            Additional keyword arguments to pass to `add_maf_file`.

        Returns
        -------
        None

        Notes
        -----
        - This function loops through the list of files and calls `add_maf_file` for each.
        - Handles missing values and ensures column consistency across files.
        """
        for f in fnames:
            self.add_maf_file(f, **kwargs)

    @require_file_exists
    def add_maf_file(
        self,
        fname,
        columns=None
    ):
        """
        Add a single MAF (Mutation Annotation Format) file to the SNV index and associated metadata, and update the index log.

        This function reads a MAF file, extracts SNV-related information,
        and adds it to the HDF5 dataset. If missing values or required columns
        are not present, they are handled accordingly.

        Parameters
        ----------
        fname : str
            Path to the MAF file to be processed.
        columns : dict, optional
            Mapping of file column names to SNV metadata columns.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If an error occurs during processing, the HDF5 file is closed, and an exception is raised.

        Notes
        -----
        - Column names expected are GDC MAF Format v1.0.0. If provided MAF files have a different format,
          the `columns` argument should be used to map column names to DNAStream SNV metadata columns.
        - SNV labels are created by concatenating the first four columns (chr:pos:ref:alt).
        - The extracted data is stored in the HDF5 file in a structured format.
        """
        if columns:
            column_dict = columns 
        else:
            column_dict = {
                "Chromosome": "chrom",
                "Start_Position": "pos",
                "End_Position": "end_pos",
                "Reference_Allele": "ref_allele",
                "Tumor_Seq_Allele2": "alt_allele",
                "Hugo_Symbol": "hugo",
                "Entrez_Gene_Id": "gene",
                "Filter" : "filter",
                "Variant_Classification": "variant_classification",
                "Variant_Type" : "variant_type",
                "dbSNP_RS" : "dbsnp_id"
            }

    
        try:
            self.load_metadata(fname, index_name=GlobalIndexName.SNV.value,delimiter="\t",
                                label_col=["chrom", "pos", "ref_allele", "alt_allele"], 
                                label_sep=":", columns = column_dict)

        except Exception:
            raise
    
    @require_file_exists
    def load_metadata(self, fname, index_name, label_col, delimiter=",", label_sep=":", columns=None):
        """
        Read sample metadata from a CSV file, add any new sample names to the index,
        and insert metadata into the /sample/metadata table in the HDF5 file.

        Parameters
        ----------
        fname : str
            Path to the metadata CSV file.
        index_name : str
            Name of the index to be updated.
        label_col : str or list of str
            Column name(s) in the CSV that contain the labels for the index.
        delimiter : str, optional
            Delimiter used in the metadata input file (default is ",").
        label_sep : str, optional
            Separator used to concatenate label columns if label_col is a list (default is ":").
        columns : dict, optional
            Mapping of column names in file to column names in metadata.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the index_name is not found in global_idx.
        Exception
            For other errors during file processing.
        """
        try:
            metadata_dict = {}
           
            with open(fname, newline="") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                if not isinstance(label_col, list):
                    label_col = [label_col]
                if columns is not None:
                    reader.fieldnames = [columns.get(col, col) for col in reader.fieldnames]
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
            raise

    # -----------------------------------
    # Copy Number I/O Functions
    # -----------------------------------
    def add_ascat_sc_copy_numbers_from_directory(
        self,
        directory,
        file_regex="*.ASCAT.scprofile.txt",
        sample_label_regex=r".*(?=\.bam)",
        allele_specific=True,
        modality=Modalities.SCDNA.value,
    ):
        """
        Load ASCAT SC copy number files from a directory and extract sample labels from filenames.

        Parameters
        ----------
        directory : str
            Directory containing ASCAT SC files.
        file_regex : str, optional
            Glob-style pattern to select files (default: '*.ASCAT.scprofile.txt').
        sample_label_regex : str, optional
            Regular expression to extract sample label from filename. The entire match (group 0) will be used as the label.
        allele_specific : bool, optional
            Whether to load allele-specific copy numbers (default: True).
        modality : str, optional
            Modality to assign to each sample (default: 'scdna').

        Returns
        -------
        list of str
            List of file paths that were skipped due to errors or missing sample label matches.

        Raises
        ------
        Exception
            If there is an error loading a file, it is added to the skipped_files list.
        """

        if allele_specific:
            ascat_sc_fn = self.add_ascat_sc_allele_specific_copy_numbers
        else:
            ascat_sc_fn = self.add_ascat_sc_total_copy_numbers

        skipped_files = []
        pattern = os.path.join(directory, file_regex)
        file_list = glob.glob(pattern)

        for f in file_list:
            
            basename = os.path.basename(f)
            match = re.search(sample_label_regex, basename)
            if not match:
                if self.verbose:
                    print(f"#Skipping {basename}: no match for regex '{sample_label_regex}'")
                continue

            sample_label = match.group(1) if match.lastindex else match.group(0)
            # if self.verbose:
            
            #     print(f"#Loading ASCAT SC file: {basename} → sample_label: {sample_label}")
            try:
                ascat_sc_fn(f, sample_label, modality=modality)
            except Exception as e:
                skipped_files.append(f)
                if self.verbose:
                    print(f"#Failed to load {f} due to error: {e}")
        return skipped_files


    def add_ascat_sc_total_copy_numbers(
        self, fname, sample_label, modality=Modalities.SCDNA.value, columns=None
    ):
        """
        Parse an ASCAT SC total copy number file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the ASCAT SC total copy number file to parse.
        sample_label : str
            Label for the sample being added.
        modality : str, optional
            Modality of the sample (e.g., "scdna", "lcm"). Default is "scdna".
        columns : dict, optional
            Mapping of file column names to required column names (["chromosome", "start", "end", "logr", "total_copy_number"]).

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the required columns are missing or modality is not recognized.
        Exception
            If file parsing or data insertion fails.
        """

        try:

            try:
                modality_enum = Modalities(modality)
            except ValueError:
                raise ValueError(f"Modality '{modality}' not recognized.")

            df = pd.read_csv(fname, sep="\t")
        
            if columns is not None:
                df= df.rename(columns=columns)
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
                _ = self.copy_number_scdna_view.register([sample_label])
            elif modality_enum == Modalities.LCM:
                _ = self.copy_number_lcm_view.register([sample_label])


            sample_idx = self.copy_number_scdna_view.label_to_view_idx(sample_label)


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
            raise

    @require_file_exists
    def add_ascat_sc_allele_specific_copy_numbers(
        self, fname, sample_label, modality=Modalities.SCDNA.value, columns=None
    ):
        """
        Parse an ASCAT SC allele-specific copy number file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the ASCAT SC allele-specific copy number file to parse.
        sample_label : str
            Label for the sample being added.
        modality : str, optional
            Modality of the sample (e.g., "scdna", "lcm"). Default is "scdna".
        columns : dict, optional
            Mapping of file column names to required column names (["chromosome", "start", "end", "logr", "total_copy_number"]).



        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the required columns are missing or modality is not recognized.
        Exception
            If file parsing or data insertion fails.
        """
    
        try:
    
            try:
                modality_enum = Modalities(modality)
            except ValueError:
                raise ValueError(f"Modality '{modality}' not recognized.")

            df = pd.read_csv(fname, sep="\t")

            if columns is not None:
                df= df.rename(columns=columns)
            

            required_columns = ["chr", "startpos", "endpos", "logr", "BAF", "nA", "nB"]
            missing = [col for col in required_columns if col not in df.columns]
            if missing:
                raise ValueError(f"Missing columns in {fname}: {missing}")

            self.batch_sample_add([sample_label], source_file=fname)
            self.insert_sample_metadata(
                {sample_label: {"modality": modality_enum.value}}
            )

            if modality_enum == Modalities.SCDNA:
                label_to_idx = self.copy_number_scdna_view.register([sample_label])
            elif modality_enum == Modalities.LCM:
                label_to_idx =self.copy_number_lcm_view.register([sample_label])

            sample_idx = label_to_idx[sample_label]
 
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
            print("addding data")
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
            raise


    @require_file_exists
    def add_battenberg_copy_numbers(self, fname, sample_label):
        """
        Parse a Battenberg file and extract the copy number data.

        Parameters
        ----------
        fname : str
            Path to the Battenberg file to parse.
        sample_label : str
            The sample label of the file.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If file parsing or data insertion fails.
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
            
           

            label_to_idx=  self.copy_number_bulk_view.register([sample_label])

            sample_idx = self.copy_number_bulk_view.label_to_view_idx(sample_label)
           
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
            raise

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
                f"Shape of `data` {dat.shape} does not match seg_indices {len(seg_indices)} x sample_indices {len(sample_indices)}"
            )

     
        for i, seg in enumerate(seg_indices):
            for j, sample in enumerate(sample_indices):
                try:
                    entry = dat[i, j]
                    if entry.dtype != ALLELE_SPECIFIC_CN_DTYPE:
                        entry = np.asarray(entry, dtype=ALLELE_SPECIFIC_CN_DTYPE)
                    self[table_name][seg, sample] = entry
                except Exception as e:
                    print(f"Failed to write to ({seg}, {sample}) with entry: {entry}")
                    print(f"Entry type: {type(entry)}, dtype: {getattr(entry, 'dtype', 'N/A')}")
                    raise e

    # -----------------------------------
    # Tree I/O functions
    # -----------------------------------

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


    @require_file_exists
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

    # -----------------------------------
    # Read Count I/O functions
    # -----------------------------------
    @require_file_exists
    def add_read_counts(self, fname, columns=None):
        """
        Add variant and total read count data for SNVs from a file.

        This method reads a file containing variant and total read counts, adds SNVs and samples to the index,
        maps SNV and sample labels to indices, and updates the HDF5 dataset efficiently.

        Parameters
        ----------
        fname : str
            Path to the input file containing read counts.
        columns : dict, optional
            Mapping of the column names in the file to the expected columns ["snv", "sample", "alt", "total"].

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If required columns are missing.
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

       
            for arr in ["variant", "total"]:
                self._log_dataset_modification(
                    f"read_counts/{arr}", operation="update", source_file=fname
                )

        except Exception:
            raise
