import pytest
from dnastream.provenance import Provenance
from dnastream._builtin_schemas import PROVENANCE_LOG_SCHEMA


"""
-Test that exception is raised if provenance tries to log itself (done)
-Test that events are added to the log properly via DNAStream (done)
-Test that exception is raised if hooks are registered on Provenance objects (done)
"""


def test_provenance_correct_group_enforcement(temp_h5_handle):
    """
    Given an existing handle to a group
    when the group name is not "/registry" initialization
    of the registry is prohibited and a ValueError is thrown,
    when the group name is "/registry"
    """
    grp = temp_h5_handle.require_group("foo")
    with pytest.raises(ValueError):
        Provenance(grp, "test")

    grp = temp_h5_handle.require_group("provenance")
    Provenance(grp, "test")


def test_dataset_not_created(temp_h5_handle):
    grp = temp_h5_handle.require_group("provenance")
    log = Provenance(grp, "test")
    with pytest.raises(RuntimeError):
        log.validate(strict=True)


def test_provenance_does_not_log_itself(temp_h5_handle):
    grp = temp_h5_handle.require_group("provenance")
    log = Provenance(grp, "test")
    log.create(schema=PROVENANCE_LOG_SCHEMA)
    assert len(log) == 0

    with pytest.raises(RuntimeError):
        log.register_hook()
        log._emit()


def test_provenance_is_valid(temp_h5_handle):
    grp = temp_h5_handle.require_group("provenance")
    log = Provenance(grp, "test")
    log.create(PROVENANCE_LOG_SCHEMA)
    log.validate()


def test_proper_logging(dnastream_obj, temp_sample_csv):
    assert len(dnastream_obj.log) > 0
    init_len = len(dnastream_obj.log)
    dnastream_obj.io.add_samples_from_files(temp_sample_csv)
    assert len(dnastream_obj.log) == init_len + 2

    log_df = dnastream_obj.log.to_dataframe()
    # print(log_df.event)

    assert log_df.shape[0] == len(dnastream_obj.log)
