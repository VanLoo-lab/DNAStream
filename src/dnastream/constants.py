import h5py
from typing import Literal, Callable
import numpy as np
import uuid
from enum import Enum


SCHEMA_VERSION: str = "0.1.1"


# Provenance constants
# EVENTS: tuple[str, ...] = ("append", "create", "modify", "designate")
class EVENTS(str, Enum):
    APPEND = "append"
    CREATE = "create"
    MODIFY = "modify"
    DESIGNATE = "designate"


class SCOPE(str, Enum):
    DNASTREAM = "DNAStream"
    REGISTRY = ("registry",)
    IO = "io"


# Schema Constants
STR_DTYPE = h5py.string_dtype("utf-8")

# Registry constants
VALID_MODES = {"active_only", "all", "non_active"}  # registry modes
BY_OPTIONS = {"label", "id"}
REGISTRY_SPINE = (
    "active",
    "id",
    "idx",
    "label",
    "created_at",
    "created_by",
    "modified_by",
    "modified_at",
)

LABEL_SCALARS = (str, bytes, np.str_, np.bytes_)  # allow only label-ish things
ID_SCALARS = (str, bytes, np.str_, np.bytes_, uuid.UUID)


# typing
Event = Literal[EVENTS.APPEND, EVENTS.CREATE, EVENTS.DESIGNATE, EVENTS.MODIFY]
Scope = Literal[SCOPE.DNASTREAM, SCOPE.REGISTRY, SCOPE.IO]
Hook = Callable[..., None]
