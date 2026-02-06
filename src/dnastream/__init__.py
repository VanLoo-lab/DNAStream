from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dnastream")  # must match pyproject.toml [project].name
except PackageNotFoundError:
    # Running from a source tree without an installed distribution
    __version__ = "0+unknown"

from .dnastream import DNAStream

__all__ = ["DNAStream", "__version__"]
