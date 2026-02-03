import pytest
from dnastream.registry import Registry
import numpy as np
import pandas
import uuid


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


def test_fields_property(registry_obj):
    expected = {"variable", "value"}
    assert set(registry_obj.fields) == expected


def test_registry_not_created_twice(registry_obj, temp_registry_schema):
    """
    Test to make sure the registry is not recreated if it already exists
    """
    with pytest.raises(RuntimeError):
        registry_obj.create(temp_registry_schema, if_exists="raise")

    with pytest.raises(ValueError):
        registry_obj.create(temp_registry_schema, if_exists="foo")

    with pytest.warns(UserWarning):
        registry_obj.create(temp_registry_schema, if_exists="open")


# # test that the validate function catches
# """
#   - checks that ids are unique
#     - checks that registry path is valid
#     - checks that fields include required fields id, label, idx and active
#     - checks row idx are valid, i.e., between 0...size -1 and idx equals row index
#     - checks ids are non-null and valid uuids
#     - checks labels are unique in active only
# """


def test_registry_is_open_with_cache_mode(
    registry_obj,
):
    """
    - test that rows are added to the registry
    -test that new unique labels are activated
    - test that collisions respect the value of keep_new
    - test that the data size on disk matches the size of the input
    - test that an entry on disk matches the input

    """
    _ = registry_obj.open()
    assert registry_obj._cache_valid
    assert isinstance(registry_obj._label_to_idx, dict)
    assert isinstance(registry_obj._id_to_idx, dict)


def test_add_unique_rows_to_registry(
    registry_obj, temp_data_rows, temp_registry_schema
):
    """
    - test that rows are added to the registry on disk
     -test that new unique labels are activated
    - test that the data size on disk matches the size of the input
    - test that an entry on disk matches the input

    """
    registry_obj.add(temp_data_rows)
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]
    assert len(registry_obj) == len(temp_data_rows)
    assert not registry_obj._cache_valid

    registry_obj.validate()
    ids = registry_obj.resolve_ids(labs)
    print(ids)

    assert len(ids) == len(temp_data_rows)

    disk_labels = registry_obj.resolve_labels(ids)
    print(disk_labels)

    assert len(ids) == len(labs)
    assert len(disk_labels) == len(temp_data_rows)
    assert set(disk_labels) == set(labs)

    # check that all the labels have been activated
    with pytest.warns(UserWarning):
        df = registry_obj.get(labs, by="label", mode="non_active")
    assert df.shape[0] == 0

    df = registry_obj.get(labs, by="label", mode="all")
    assert np.all(df["active"])
    assert df.shape[0] == len(registry_obj)

    # # should return only 1 row with the "value" stored as 0 on disk
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
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]
    registry_obj.add(temp_data_rows)
    assert len(registry_obj) == len(temp_data_rows)
    assert not registry_obj._cache_valid

    df = registry_obj.get(labs, by="label", mode="active_only")
    assert df.shape[0] == len(temp_data_rows)

    registry_obj.validate()

    registry_obj.add(temp_data_rows, activate_new=False, allow_duplicate_labels=True)
    assert len(registry_obj) == 2 * len(temp_data_rows)
    assert not registry_obj._cache_valid

    registry_obj.validate()

    registry_obj.add(temp_data_rows, activate_new=True, allow_duplicate_labels=True)
    assert len(registry_obj) == 3 * len(temp_data_rows)
    df = registry_obj.get(labs, by="label", mode="active_only")
    registry_obj.validate()

    df = registry_obj.get(labs, by="label", mode="active_only")
    assert df.shape[0] == len(temp_data_rows)


def test_find_ids(registry_obj, temp_data_rows, temp_registry_schema):
    """
    Test that find_ids correctly returns multiple ids when there are duplicate labels
    """
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    row_vals = registry_obj.find_ids(labs, mode="all")
    for _, rows in row_vals.items():
        assert len(rows) == 1

    registry_obj.add(temp_data_rows, allow_duplicate_labels=True)
    row_vals = registry_obj.find_ids(labs, mode="all")
    for _, rows in row_vals.items():
        assert len(rows) == 2


def test_to_dataframe(registry_obj, temp_data_rows):
    """
    Docstring for test_to_dataframe

    :param registry_obj: Description
    :param temp_data_rows: Description
    """
    registry_obj.add(temp_data_rows)

    df = registry_obj.to_dataframe(mode="all")
    assert isinstance(df, pandas.DataFrame)
    assert df.shape[0] == len(temp_data_rows)

    registry_obj.add(temp_data_rows, allow_duplicate_labels=True)
    df = registry_obj.to_dataframe(mode="active_only")
    assert df.shape[0] == len(temp_data_rows)

    df = registry_obj.to_dataframe(mode="all")
    assert df.shape[0] == 2 * len(temp_data_rows)

    variables = df["variable"].tolist()
    # variables = [s.decode("utf-8") for s in variables]
    values = df["value"].to_list()

    for var, val in zip(variables, values):
        assert int(var.replace("var", "")) == val


def test_contains(registry_obj, temp_data_rows):
    """
    Test contains as a dunder method
    """
    registry_obj.add(temp_data_rows)
    df = registry_obj.to_dataframe()
    ids = df["id"].tolist()
    for i in ids:
        assert i in registry_obj

    assert uuid.uuid4() not in registry_obj


def test_get_item(registry_obj, temp_data_rows):
    """
    Test __getitem__ dunder method
    """
    with pytest.raises(KeyError):
        bad_id = uuid.uuid4()
        _ = registry_obj[bad_id]

    registry_obj.add(temp_data_rows)


#     df = registry_obj.to_dataframe()
#     ids = df["id"].tolist()
#     for i in ids:
#         entity = registry_obj[i]
#         variable = entity["variable"]
#         intvar = int(variable.replace("var", ""))
#         assert  intvar == entity["value"]


def test_iter(registry_obj, temp_data_rows):
    """
    Test iter dunder method
    """
    registry_obj.add(temp_data_rows)
    for i, row in enumerate(registry_obj):
        print(row["variable"])

        assert row["value"] == temp_data_rows[i]["value"]
        assert row["variable"] == temp_data_rows[i]["variable"]


def test_resolve_ids(registry_obj, temp_data_rows, temp_registry_schema):
    """
    Test contains dunder method
    """
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]
    registry_obj.add(temp_data_rows)
    reg_id = registry_obj.resolve_ids(labs[0])
    assert reg_id is not None
    reg_ids = registry_obj.resolve_ids(labs)
    ids = []
    for i, row in enumerate(registry_obj):
        reg_ids[i] = row["id"]
        ids.append(row["id"])

    # test to make sure missing labels return missing field
    id = registry_obj.resolve_ids("foo", missing=None)

    assert id is None

    # #ensure lookup only returns the active label
    # registry_obj.add(temp_data_rows, temp_registry_schema, activate_new=False)
    # ids = registry_obj.resolve_ids(labs)
    # assert len(ids) ==len(temp_data_rows)
    # assert ids[0] == id


def test_resolve_labels(registry_obj, temp_data_rows, temp_registry_schema):
    """
    Test contains as a dunder method
    """
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    ids = [row["id"] for row in registry_obj]

    reg_labels = registry_obj.resolve_labels(ids)
    for reg_label, lab in zip(reg_labels, labs):
        assert reg_label == lab

    # test to make sure missing labels return missing field
    reg_label = registry_obj.resolve_labels(uuid.uuid4(), missing=None)
    assert reg_label is None

    reg_label = registry_obj.resolve_labels(str(uuid.uuid4()), missing=None)
    assert reg_label is None

    # ensure resolve only returns the active label
    registry_obj.add(temp_data_rows, activate_new=False)
    reg_labels = registry_obj.resolve_labels(ids)
    assert len(labs) == len(reg_labels)
    assert labs[0] == reg_labels[0]


def test_activate_ids(registry_obj, temp_data_rows, temp_registry_schema):

    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    first_ids = registry_obj.resolve_ids(labs)
    registry_obj.add(temp_data_rows)
    ids = [row["id"] for row in registry_obj]
    ids = ids[::-1]
    activate_ids = ids[: len(temp_data_rows)]
    registry_obj.activate_ids(activate_ids)
    df = registry_obj.to_dataframe()

    # validate that newly added ids with duplicate labels are active
    all_active = df[df["id"].isin(activate_ids)]
    assert np.all(all_active["active"])

    assert np.sum(df["active"]) == len(temp_data_rows)

    # make sure duplicate labels are deactivated
    not_active = df[df["id"].isin(first_ids)]
    assert np.all(~not_active["active"])
    registry_obj.validate()


def test_deactivate_ids(registry_obj, temp_data_rows, temp_registry_schema):
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    first_ids = registry_obj.resolve_ids(labs)
    registry_obj.deactivate_ids(first_ids)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == 0

    registry_obj.add(temp_data_rows)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == len(temp_data_rows)

    registry_obj.deactivate_ids(first_ids)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == len(temp_data_rows)
    registry_obj.validate()


def test_deactivate_ids(registry_obj, temp_data_rows, temp_registry_schema):
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    # add all temp data and deactivate ids them
    registry_obj.add(temp_data_rows)
    first_ids = registry_obj.resolve_ids(labs)
    registry_obj.deactivate_ids(first_ids)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == 0

    # add the duplicate data and ensure the default activate policy works
    registry_obj.add(temp_data_rows)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == len(temp_data_rows)

    with pytest.warns(UserWarning, match="not found|Requested selection"):
        registry_obj.deactivate_ids(labs, warn_missing=True)

    with pytest.raises(ValueError):

        registry_obj.deactivate_ids([i for i in range(5)])


def test_activate_ids(registry_obj, temp_data_rows, temp_registry_schema):

    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    first_ids = registry_obj.resolve_ids(labs)
    registry_obj.add(temp_data_rows, allow_duplicate_labels=True)
    ids = [row["id"] for row in registry_obj]
    ids = ids[::-1]
    activate_ids = ids[: len(temp_data_rows)]
    registry_obj.activate_ids(activate_ids)
    df = registry_obj.to_dataframe()

    # validate that newly added ids with duplicate labels are active
    all_active = df[df["id"].isin(activate_ids)]
    assert np.all(all_active["active"])

    assert np.sum(df["active"]) == len(temp_data_rows)

    # make sure duplicate labels are deactivated
    not_active = df[df["id"].isin(first_ids)]
    assert np.all(~not_active["active"])
    registry_obj.validate()


def test_deactivate_labels(registry_obj, temp_data_rows, temp_registry_schema):
    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    registry_obj.add(temp_data_rows)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == len(temp_data_rows)

    # check that value error when invalid label types are passed
    with pytest.raises(ValueError):
        registry_obj.deactivate_labels([i for i in range(5)])

    with pytest.warns(UserWarning):
        registry_obj.deactivate_ids("var10")

    registry_obj.deactivate_labels(labs)
    df = registry_obj.to_dataframe()
    assert np.sum(df["active"]) == 0
    registry_obj.validate()


def test_activate_labels(registry_obj, temp_data_rows, temp_registry_schema):

    labs = [temp_registry_schema.make_label(mydict) for mydict in temp_data_rows]

    with pytest.warns(UserWarning):
        registry_obj.activate_labels(["var10"], warn_missing=True)

    registry_obj.add(temp_data_rows)
    df = registry_obj.to_dataframe()
    assert df["active"].sum() == len(temp_data_rows)

    registry_obj.deactivate_labels(labs)
    df = registry_obj.to_dataframe()
    assert df["active"].sum() == 0

    registry_obj.activate_labels(labs)
    df = registry_obj.to_dataframe()
    assert df["active"].sum() == len(temp_data_rows)
    registry_obj.add(temp_data_rows, allow_duplicate_labels=True)

    # check that the newest labels are activated
    registry_obj.activate_labels(labs, activate_newest=True)
    registry_obj.activate_labels(labs)
    df = registry_obj.to_dataframe()
    assert df["active"].sum() == len(temp_data_rows)
    active_states = df["active"].to_list()
    assert sum([active_states[i] for i in range(len(temp_data_rows))]) == 0
    assert sum(
        [active_states[len(temp_data_rows) + i] for i in range(len(temp_data_rows))]
    ) == len(temp_data_rows)

    with pytest.raises(ValueError):
        registry_obj.activate_labels([i for i in range(5)])

    registry_obj.validate()


def test_update_registry_items(registry_obj, temp_data_rows):
    """Test editable and only editable metadata data fields are updated."""
    registry_obj.add(temp_data_rows)

    df = registry_obj.to_dataframe()
    ids = df["id"].to_numpy()
    # ids = registry_obj.resolve_labels(labs)

    updates = [{"id": i, "value": 5} for i in ids]

    registry_obj.update(updates)
    df = registry_obj.to_dataframe()
    values = df["value"].to_numpy()
    assert np.all(values == 5)

    # test the size of the registy  hasn't changed
    assert len(registry_obj) == len(temp_data_rows)

    key = uuid.uuid4()
    print(f"fake key: {key}")
    bad_id_update = [{"id": key, "value": 5}]
    non_editable_field_update = [{"id": ids[0], "label": "foo", "value": 10}]
    label_from_field_update = [{"id": ids[0], "variable": "foo"}]
    with pytest.warns(UserWarning):
        registry_obj.update(bad_id_update, warn_missing=True)
        registry_obj.update(non_editable_field_update, warn_missing=False)
        registry_obj.update(label_from_field_update, warn_missing=False)

    # check that the value was updated but the label wasn't
    row_dict = registry_obj[ids[0]]
    assert row_dict["value"] == 10
    # test
    missing_update = [{"value": 5}]
    with pytest.raises(ValueError):
        registry_obj.update(missing_update, warn_missing=True)

    # with pytest.raises(KeyError):
    #     registry_obj.update(missing_update, temp_registry_schema, warn_missing=True)

    # registry_obj.update(updates, temp_registry_schema)

    # assert len(registry_obj) == len(temp_data_rows)
    # assert not registry_obj._cache_valid

    # registry_obj.validate()
    # ids = registry_obj.resolve_ids(labs)
    # assert len(ids) == len(temp_data_rows)
    # disk_labels = registry_obj.resolve_labels(ids)
    # assert len(ids) == len(labs)
    # assert len(disk_labels) == len(temp_data_rows)
    # assert set(disk_labels) == set(labs)

    # # check that all the labels have been activated
    # df = registry_obj.get(labs, by="label", mode="non_active")
    # assert df.shape[0] == 0

    # df = registry_obj.get(labs, by="label", mode="all")
    # assert np.all(df["active"])
    # assert df.shape[0] == len(registry_obj)

    # # should return only 1 row with the "value" stored as 0 on disk
    # df = registry_obj.get("var0", by="label")
    # assert df.shape[0] == 1
    # assert df["value"][0] == 0


def test_add_duplicate_labels_false(registry_obj, temp_data_rows):
    """
    Test to make sure that setting `allow_duplicate_labels` to false does not
    add duplicate labels to the registry.
    """

    registry_obj.add(temp_data_rows, activate_new=True)
    assert len(registry_obj) == len(temp_data_rows)

    registry_obj.add(temp_data_rows, activate_new=True, allow_duplicate_labels=False)
    assert len(registry_obj) == len(temp_data_rows)

    registry_obj.add(temp_data_rows, activate_new=False, allow_duplicate_labels=False)
    assert len(registry_obj) == len(temp_data_rows)
