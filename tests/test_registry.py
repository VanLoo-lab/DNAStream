import pytest
from dnastream.registry import Registry
import numpy as np
import pandas


def test_registry_correct_group_enforcement(temp_h5_handle):
    """
    Given an existing handle to a group
    when the group name is not "/registry" initialization
    of the registry is prohibited and a ValueError is thrown,
    when the group name is "/registry"
    """
    grp = temp_h5_handle.require_group("foo")
    with pytest.raises(ValueError):
        Registry(grp, "test")

    grp = temp_h5_handle.require_group("registry")
    Registry(grp, "test")


# test that the validate function catches
"""
  - checks that ids are unique
    - checks that registry path is valid
    - checks that fields include required fields id, label, idx and active
    - checks row idx are valid, i.e., between 0...size -1 and idx equals row index
    - checks ids are non-null and valid uuids
    - checks labels are unique in active only
"""


def test_registry_is_open_with_cache_mode(
    registry_obj, temp_registry_schema, temp_data_rows
):
    """
    - test that rows are added to the registry
    -test that new unique labels are activated
    - test that collisions respect the value of keep_new
    - test that the data size on disk matches the size of the input
    - test that an entry on disk matches the input

    """
    _ = registry_obj.open(temp_registry_schema, mode="active_only")
    assert registry_obj._cache_valid
    assert registry_obj._cache_mode == "active_only"
    assert isinstance(registry_obj._label_to_idx, dict)
    assert isinstance(registry_obj._id_to_idx, dict)

    _ = registry_obj.open(temp_registry_schema, mode="all")
    assert registry_obj._cache_valid
    assert registry_obj._cache_mode == "all"
    assert registry_obj._label_to_idx is None
    assert isinstance(registry_obj._id_to_idx, dict)

    # test won't open in invalid mode

    with pytest.raises(ValueError):
        _ = registry_obj.open(temp_registry_schema, mode="foo")


def test_add_unique_rows_to_registry(
    registry_obj, temp_data_rows, temp_registry_schema
):
    """
    - test that rows are added to the registry on disk
     -test that new unique labels are activated
    - test that the data size on disk matches the size of the input
    - test that an entry on disk matches the input

    """
    registry_obj.add(temp_data_rows, temp_registry_schema)
    labs = [
        temp_registry_schema["label_normalizer"](mydict["variable"])
        for mydict in temp_data_rows
    ]
    assert len(registry_obj) == len(temp_data_rows)
    assert not registry_obj._cache_valid

    registry_obj.validate()
    ids = registry_obj.lookup_id(labs)
    assert len(ids) == len(temp_data_rows)
    disk_labels = registry_obj.lookup_label(ids)
    assert len(ids) == len(labs)
    assert len(disk_labels) == len(temp_data_rows)
    assert set(disk_labels) == set(labs)

    # check that all the labels have been activated
    df = registry_obj.get(labs, by="label", mode="non_active")
    assert df.shape[0] == 0

    df = registry_obj.get(labs, by="label", mode="all")
    assert np.all(df["active"])
    assert df.shape[0] == len(registry_obj)

    # should return only 1 row with the "value" stored as 0 on disk
    df = registry_obj.get("var0", by="label")
    assert df.shape[0] == 1
    assert df["value"][0] == 0


def test_add_duplicate_rows_activated_properly(
    registry_obj, temp_data_rows, temp_registry_schema
):
    """
    Test that new rows with duplicate labels are actived/activated properly according to the Registry contract.
    - test that new rows are activated by default
    - test that old rows with dup labels stay activated and new rows are added but not activated when activate_new=False
    - test that new rows are activated when dup labels exist and old rows with duplicate labels are deactivated
    """
    labs = [
        temp_registry_schema["label_normalizer"](mydict["variable"])
        for mydict in temp_data_rows
    ]

    registry_obj.add(temp_data_rows, temp_registry_schema)
    assert len(registry_obj) == len(temp_data_rows)
    assert not registry_obj._cache_valid
    df = registry_obj.get(labs, by="label", mode="active_only")
    assert df.shape[0] == len(temp_data_rows)

    registry_obj.validate()

    registry_obj.add(temp_data_rows, temp_registry_schema, activate_new=False)
    assert len(registry_obj) == 2 * len(temp_data_rows)
    assert not registry_obj._cache_valid

    registry_obj.validate()

    registry_obj.add(temp_data_rows, temp_registry_schema, activate_new=True)
    assert len(registry_obj) == 3 * len(temp_data_rows)
    df = registry_obj.get(labs, by="label", mode="active_only")
    registry_obj.validate()

    df = registry_obj.get(labs, by="label", mode="active_only")
    assert df.shape[0] == len(temp_data_rows)


def test_find_ids(registry_obj, temp_data_rows, temp_registry_schema):
    """
    Test that find_ids correctly returns multiple ids when there are duplicate labels
    """
    labs = [
        temp_registry_schema["label_normalizer"](mydict["variable"])
        for mydict in temp_data_rows
    ]

    registry_obj.add(temp_data_rows, temp_registry_schema)
    row_vals = registry_obj.find_ids(labs, mode="all")
    for _, rows in row_vals.items():
        assert len(rows) == 1

    registry_obj.add(temp_data_rows, temp_registry_schema)
    row_vals = registry_obj.find_ids(labs, mode="all")
    for _, rows in row_vals.items():
        assert len(rows) == 2


def test_to_dataframe(registry_obj, temp_data_rows, temp_registry_schema):
    """
    Docstring for test_to_dataframe

    :param registry_obj: Description
    :param temp_data_rows: Description
    :param temp_registry_schema: Description
    """
    registry_obj.add(temp_data_rows, temp_registry_schema)

    df = registry_obj.to_dataframe(mode="all")
    assert isinstance(df, pandas.DataFrame)
    assert df.shape[0] == len(temp_data_rows)

    registry_obj.add(temp_data_rows, temp_registry_schema)
    df = registry_obj.to_dataframe(mode="active_only")
    assert df.shape[0] == len(temp_data_rows)

    df = registry_obj.to_dataframe(mode="all")
    assert df.shape[0] == 2 * len(temp_data_rows)

    variables = df["variable"].tolist()
    variables = [s.decode("utf-8") for s in variables]
    values = df["value"].to_list()

    for var, val in zip(variables, values):
        assert int(var.replace("var", "")) == val


# def test_lookup_labels(registry_obj, temp_data_rows, temp_registry_schema):
#     """
#     Test that find_ids correctly returns multiple ids when there are duplicate labels
#     """
#     labs = [temp_registry_schema["label_normalizer"](mydict["variable"]) for mydict in temp_data_rows]

#     registry_obj.add(temp_data_rows, temp_registry_schema)
#     df = registry_obj.to_dataframe()

#     registry_obj.add(temp_data_rows, temp_registry_schema)
#     row_vals = registry_obj.find_ids(labs, mode="all")
#     for _, rows in row_vals.items():
#         assert len(rows) ==2
