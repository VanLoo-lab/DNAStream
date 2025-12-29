from dnastream._h5base import H5Dataset


class _TestH5Dataset(H5Dataset):
    """
    Helper class used to test the ABC H5Dataset
    """

    def add(self, *args, **kwargs):
        raise NotImplementedError
