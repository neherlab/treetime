from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class node(_message.Message):
    __slots__ = ("mutation_positions", "mutation_other_fields", "ignored_range_start", "ignored_range_end", "node_id", "children_offsets", "children_lengths", "condensed_nodes", "changed")
    MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    IGNORED_RANGE_START_FIELD_NUMBER: _ClassVar[int]
    IGNORED_RANGE_END_FIELD_NUMBER: _ClassVar[int]
    NODE_ID_FIELD_NUMBER: _ClassVar[int]
    CHILDREN_OFFSETS_FIELD_NUMBER: _ClassVar[int]
    CHILDREN_LENGTHS_FIELD_NUMBER: _ClassVar[int]
    CONDENSED_NODES_FIELD_NUMBER: _ClassVar[int]
    CHANGED_FIELD_NUMBER: _ClassVar[int]
    mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    ignored_range_start: _containers.RepeatedScalarFieldContainer[int]
    ignored_range_end: _containers.RepeatedScalarFieldContainer[int]
    node_id: int
    children_offsets: _containers.RepeatedScalarFieldContainer[int]
    children_lengths: _containers.RepeatedScalarFieldContainer[int]
    condensed_nodes: _containers.RepeatedScalarFieldContainer[str]
    changed: int
    def __init__(self, mutation_positions: _Optional[_Iterable[int]] = ..., mutation_other_fields: _Optional[_Iterable[int]] = ..., ignored_range_start: _Optional[_Iterable[int]] = ..., ignored_range_end: _Optional[_Iterable[int]] = ..., node_id: _Optional[int] = ..., children_offsets: _Optional[_Iterable[int]] = ..., children_lengths: _Optional[_Iterable[int]] = ..., condensed_nodes: _Optional[_Iterable[str]] = ..., changed: _Optional[int] = ...) -> None: ...

class node_idx(_message.Message):
    __slots__ = ("node_id", "node_name")
    NODE_ID_FIELD_NUMBER: _ClassVar[int]
    NODE_NAME_FIELD_NUMBER: _ClassVar[int]
    node_id: int
    node_name: str
    def __init__(self, node_id: _Optional[int] = ..., node_name: _Optional[str] = ...) -> None: ...

class meta(_message.Message):
    __slots__ = ("ref_nuc", "nodes_idx_next", "chromosomes", "root_offset", "root_length", "node_idx_map")
    REF_NUC_FIELD_NUMBER: _ClassVar[int]
    NODES_IDX_NEXT_FIELD_NUMBER: _ClassVar[int]
    CHROMOSOMES_FIELD_NUMBER: _ClassVar[int]
    ROOT_OFFSET_FIELD_NUMBER: _ClassVar[int]
    ROOT_LENGTH_FIELD_NUMBER: _ClassVar[int]
    NODE_IDX_MAP_FIELD_NUMBER: _ClassVar[int]
    ref_nuc: _containers.RepeatedScalarFieldContainer[int]
    nodes_idx_next: int
    chromosomes: _containers.RepeatedScalarFieldContainer[str]
    root_offset: int
    root_length: int
    node_idx_map: _containers.RepeatedCompositeFieldContainer[node_idx]
    def __init__(self, ref_nuc: _Optional[_Iterable[int]] = ..., nodes_idx_next: _Optional[int] = ..., chromosomes: _Optional[_Iterable[str]] = ..., root_offset: _Optional[int] = ..., root_length: _Optional[int] = ..., node_idx_map: _Optional[_Iterable[_Union[node_idx, _Mapping]]] = ...) -> None: ...

class sample_to_place(_message.Message):
    __slots__ = ("sample_id", "sample_mutation_positions", "sample_mutation_other_fields")
    SAMPLE_ID_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    sample_id: int
    sample_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    sample_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, sample_id: _Optional[int] = ..., sample_mutation_positions: _Optional[_Iterable[int]] = ..., sample_mutation_other_fields: _Optional[_Iterable[int]] = ...) -> None: ...

class placed_target(_message.Message):
    __slots__ = ("target_node_id", "split_node_id", "sample_mutation_positions", "sample_mutation_other_fields", "split_mutation_positions", "split_mutation_other_fields", "shared_mutation_positions", "shared_mutation_other_fields", "sample_id")
    TARGET_NODE_ID_FIELD_NUMBER: _ClassVar[int]
    SPLIT_NODE_ID_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    SPLIT_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SPLIT_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    SHARED_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SHARED_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_ID_FIELD_NUMBER: _ClassVar[int]
    target_node_id: int
    split_node_id: int
    sample_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    sample_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    split_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    split_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    shared_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    shared_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    sample_id: int
    def __init__(self, target_node_id: _Optional[int] = ..., split_node_id: _Optional[int] = ..., sample_mutation_positions: _Optional[_Iterable[int]] = ..., sample_mutation_other_fields: _Optional[_Iterable[int]] = ..., split_mutation_positions: _Optional[_Iterable[int]] = ..., split_mutation_other_fields: _Optional[_Iterable[int]] = ..., shared_mutation_positions: _Optional[_Iterable[int]] = ..., shared_mutation_other_fields: _Optional[_Iterable[int]] = ..., sample_id: _Optional[int] = ...) -> None: ...

class target(_message.Message):
    __slots__ = ("target_node_id", "parent_node_id", "sample_mutation_positions", "sample_mutation_other_fields", "split_mutation_positions", "split_mutation_other_fields", "shared_mutation_positions", "shared_mutation_other_fields")
    TARGET_NODE_ID_FIELD_NUMBER: _ClassVar[int]
    PARENT_NODE_ID_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    SPLIT_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SPLIT_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    SHARED_MUTATION_POSITIONS_FIELD_NUMBER: _ClassVar[int]
    SHARED_MUTATION_OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    target_node_id: int
    parent_node_id: int
    sample_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    sample_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    split_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    split_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    shared_mutation_positions: _containers.RepeatedScalarFieldContainer[int]
    shared_mutation_other_fields: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, target_node_id: _Optional[int] = ..., parent_node_id: _Optional[int] = ..., sample_mutation_positions: _Optional[_Iterable[int]] = ..., sample_mutation_other_fields: _Optional[_Iterable[int]] = ..., split_mutation_positions: _Optional[_Iterable[int]] = ..., split_mutation_other_fields: _Optional[_Iterable[int]] = ..., shared_mutation_positions: _Optional[_Iterable[int]] = ..., shared_mutation_other_fields: _Optional[_Iterable[int]] = ...) -> None: ...

class search_result(_message.Message):
    __slots__ = ("sample_id", "place_targets")
    SAMPLE_ID_FIELD_NUMBER: _ClassVar[int]
    PLACE_TARGETS_FIELD_NUMBER: _ClassVar[int]
    sample_id: int
    place_targets: _containers.RepeatedCompositeFieldContainer[target]
    def __init__(self, sample_id: _Optional[int] = ..., place_targets: _Optional[_Iterable[_Union[target, _Mapping]]] = ...) -> None: ...

class mutation_at_each_pos(_message.Message):
    __slots__ = ("node_id", "mut")
    NODE_ID_FIELD_NUMBER: _ClassVar[int]
    MUT_FIELD_NUMBER: _ClassVar[int]
    node_id: _containers.RepeatedScalarFieldContainer[int]
    mut: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, node_id: _Optional[_Iterable[int]] = ..., mut: _Optional[_Iterable[int]] = ...) -> None: ...

class mutation_collection(_message.Message):
    __slots__ = ("node_idx", "positions", "other_fields")
    NODE_IDX_FIELD_NUMBER: _ClassVar[int]
    POSITIONS_FIELD_NUMBER: _ClassVar[int]
    OTHER_FIELDS_FIELD_NUMBER: _ClassVar[int]
    node_idx: int
    positions: _containers.RepeatedScalarFieldContainer[int]
    other_fields: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, node_idx: _Optional[int] = ..., positions: _Optional[_Iterable[int]] = ..., other_fields: _Optional[_Iterable[int]] = ...) -> None: ...
