import pytest
from dnastream.registry import Registry


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
