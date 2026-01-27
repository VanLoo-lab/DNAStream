"""dnastream.io

I/O helpers for ingesting common genomics file formats into a
:class:`~dnastream.dnastream.DNAStream`.

The public API is intentionally thin: this module focuses on parsing files into
Python dictionaries that match registry column names, and then delegates
validation, de-duplication, and activation policies to
:class:`~dnastream.registry.Registry`.

Notes
-----
The :class:`~dnastream.io.IO` class is designed as a light wrapper around a
connected :class:`~dnastream.dnastream.DNAStream` instance.
"""

import csv
from typing import Sequence, Mapping, Any, Optional, Literal
import numpy as np
from .registry import Registry

from .utils import wrap_list, require_file_exists_static


_MAF_COLUMN_MAPPING = {
    "Chromosome": "chrom",
    "Start_Position": "start_pos",
    "End_Position": "end_pos",
    "Reference_Allele": "ref_allele",
    "Tumor_Seq_Allele2": "alt_allele",
    "Hugo_Symbol": "hugo",
    "Entrez_Gene_Id": "entrez_gene_id",
    "Filter": "filter",
    "Variant_Classification": "variant_classification",
    "Variant_Type": "variant_type",
    "dbSNP_RS": "dbsnp_id",
}


class IO:
    """File import helpers for a :class:`~dnastream.dnastream.DNAStream`.

    This class groups together convenience methods for loading external file
    formats (e.g. MAF, CSV) and appending their contents into DNAStream
    registries.

    Parameters
    ----------
    ds : dnastream.dnastream.DNAStream
        The DNAStream object that owns the target registries.

    Attributes
    ----------
    _ds : dnastream.dnastream.DNAStream
        The wrapped DNAStream instance.

    Notes
    -----
    ``IO`` is meant to be composed into
    :class:`~dnastream.dnastream.DNAStream` (e.g. as ``ds.io``) or used as a
    standalone helper that forwards attribute access to the underlying
    DNAStream.
    """

    def __init__(self, ds):
        self._ds = ds

    @property
    def ds(self):
        return self._ds

    def add_snps_from_maf(
        self,
        fname: Sequence[str] | str,
        *,
        column_mapping: Optional[Mapping[str, str]] = None,
        activate_new: bool = True,
        delimiter: str = "\t",
        **kwargs,
    ):
        """Add one or more MAF files to the 'snp' registry.

        Parameters
        ----------
        fname : list[str] or str
            Path(s) to MAF files.
        column_mapping : Mapping[str, str] or None, optional
            Mapping from file headers to registry field names.
        activate_new : bool, optional
            Passed through to Registry.add(...).
        delimiter : str, optional
            Delimiter for the input file (default tab for MAF).
        **kwargs
            Forwarded to csv.DictReader.

        Returns
        -------
        None
        """
        self._add_from_maf(
            fname=fname,
            registry_name="snp",
            column_mapping=column_mapping,
            activate_new=activate_new,
            delimiter=delimiter,
            **kwargs,
        )

    def _add_from_maf(
        self,
        fname: Sequence[str] | str,
        registry_name: Literal["snp", "variant"] = "variant",
        *,
        column_mapping: Optional[Mapping[str, str]] = None,
        activate_new: bool = True,
        delimiter: str = "\t",
        **kwargs,
    ) -> None:

        if registry_name not in {"variant", "snp"}:
            raise ValueError(
                f"Registy name '{registry_name}' invalid, must be one of 'snp' or 'variant'."
            )

        registry = getattr(self.ds, registry_name)

        if column_mapping is None:
            column_mapping = _MAF_COLUMN_MAPPING

        IO._add_files_to_registry(
            fname,
            registry,
            column_mapping=column_mapping,
            activate_new=activate_new,
            delimiter=delimiter,
            **kwargs,
        )

    def add_variants_from_maf(
        self,
        fname: Sequence[str] | str,
        *,
        column_mapping: Optional[Mapping[str, str]] = None,
        activate_new: bool = True,
        delimiter: str = "\t",
        **kwargs,
    ) -> None:
        """Add one or more MAF files to the 'variant' registry.

        Parameters
        ----------
        fname : list[str] or str
            Path(s) to MAF files.
        column_mapping : Mapping[str, str] or None, optional
            Mapping from file headers to registry field names.
        activate_new : bool, optional
            Passed through to Registry.add(...).
        delimiter : str, optional
            Delimiter for the input file (default tab for MAF).
        **kwargs
            Forwarded to csv.DictReader.

        Returns
        -------
        None
        """

        self._add_from_maf(
            fname=fname,
            registry_name="variant",
            column_mapping=column_mapping,
            activate_new=activate_new,
            delimiter=delimiter,
            **kwargs,
        )

    def add_samples_from_files(
        self,
        fname: Sequence[str] | str,
        *,
        column_mapping: Optional[Mapping[str, str]] = None,
        activate_new: bool = True,
        delimiter: str = ",",
        **kwargs,
    ) -> None:
        """Add one or more sample files to the 'sample' registry.

        Parameters
        ----------
        fname : list[str] or str
            Path(s) to sample files.
        column_mapping : Mapping[str, str] or None, optional
            Mapping from file headers to registry field names.
        activate_new : bool, optional
            Passed through to Registry.add(...).
        delimiter : str, optional
            Delimiter for the input file (default tab for MAF).
        **kwargs
            Forwarded to csv.DictReader.

        Returns
        -------
        None
        """
        registry = getattr(self.ds, "sample")
        self._add_files_to_registry(
            fname,
            registry,
            column_mapping=column_mapping,
            activate_new=activate_new,
            delimiter=delimiter,
            **kwargs,
        )

    @staticmethod
    def _coerce(name: str, v: Any, schema):
        if v is None:
            return None
        if isinstance(v, str):
            v = v.strip()
            if v == "" or v.upper() in {"NA", "N/A", "NAN", "NULL", "."}:
                return None

        if schema is None:
            return v

        try:
            field = schema.field(name)
            dt = np.dtype(field.dtype)
        except Exception:
            return v

        if dt.kind in {"i", "u"}:
            return int(v)
        if dt.kind == "f":
            return float(v)
        if dt.kind == "b":
            s = str(v).strip().lower()
            if s in {"true", "t", "1", "yes", "y"}:
                return True
            if s in {"false", "f", "0", "no", "n"}:
                return False
            return bool(v)
        return str(v)

    @staticmethod
    @require_file_exists_static
    def _parse_file(
        fname,
        columns: Sequence[str],
        column_mapping: Mapping[str, str],
        *,
        delimiter=",",
        schema=None,
        **kwargs,
    ):
        rows = []
        colset = None if columns is None else set(columns)

        with open(fname, newline="") as f:
            reader = csv.DictReader(f, delimiter=delimiter, **kwargs)

            if reader.fieldnames is None:
                raise ValueError(
                    f"File '{fname}' does not have header column. cannot parse."
                )

            if column_mapping is not None:
                reader.fieldnames = [
                    column_mapping.get(c, c) for c in reader.fieldnames
                ]

            for raw in reader:
                rec = {}
                for k, v in raw.items():

                    if colset is None or k in colset:
                        rec[k] = IO._coerce(k, v, schema)
                rows.append(rec)

        return rows

    @staticmethod
    def _add_files_to_registry(
        fname: Sequence[str] | str,
        registry: Registry,
        column_mapping: Optional[Mapping[str, str]],
        activate_new: bool = True,
        delimiter: str = "\t",
        **kwargs,
    ):

        fnames = wrap_list(fname)

        columns = registry.fields

        rows: list[dict[str, Any]] = []
        for path in fnames:

            rows.extend(
                IO._parse_file(
                    path,
                    columns,
                    column_mapping=column_mapping,
                    delimiter=delimiter,
                    schema=getattr(registry, "schema", None),
                    **kwargs,
                )
            )

        if rows:
            registry.add(rows, activate_new=activate_new)

    @staticmethod
    def parse_csv(
        fname, *, column_mapping: Mapping[str, str] | None = None, **kwargs
    ) -> list[dict]:
        """Helper method to parse a CSV file.

        Parameters
        ----------
        fname :  str
            Path(s) to sample files.
        column_mapping : Mapping[str, str] or None, optional
            Mapping from file headers to returned dictionary key names
        **kwargs
            Forwarded to csv.DictReader.

        Returns
        -------
        list[dict] : the parsed data in the file
        """
        return IO._parse_file(
            fname, columns=None, column_mapping=column_mapping, delimiter=",", **kwargs
        )

    @staticmethod
    def parse_tsv(
        fname, *, column_mapping: Mapping[str, str] | None = None, **kwargs
    ) -> list[dict]:
        """Helper method to parse a tab delimited file.

        Parameters
        ----------
        fname : str
            Path(s) to sample files.
        column_mapping : Mapping[str, str] or None, optional
            Mapping from file headers to returned dictionary key names
        **kwargs
            Forwarded to csv.DictReader.

        Returns
        -------
        list[dict] : the parsed data in the file
        """
        return IO._parse_file(
            fname, columns=None, column_mapping=column_mapping, delimiter="\t", **kwargs
        )
