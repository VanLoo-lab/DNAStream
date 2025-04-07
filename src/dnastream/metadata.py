import numpy as np
import h5py

SAMPLE_DTYPE = np.dtype(
    [
        ("sample_name", "S100"),
        ("organism", "S40"),
        ("library_strategy", "S40"),
        ("library_source", "S40"),
        ("library_selection", "S40"),
        ("library_layout", "S10"),
        ("platform", "S40"),
        ("model", "S40"),
        ("center_name", "S60"),
        ("run", "S20"),
        ("study", "S20"),
        ("coverage", "f4"),
        ("modality", "S10"),
        ("location", "S60"),
        ("bam_file_path", h5py.special_dtype(vlen=str)),
        ("batch_id", "S40"),
        ("reference_build", "S20"),
        ("pipeline_version", "S20"),
        ("date_of_sequencing", "S20"),
    ]
)
