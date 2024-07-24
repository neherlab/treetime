from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class mut(_message.Message):
    __slots__ = ("position", "ref_nuc", "par_nuc", "mut_nuc", "chromosome")
    POSITION_FIELD_NUMBER: _ClassVar[int]
    REF_NUC_FIELD_NUMBER: _ClassVar[int]
    PAR_NUC_FIELD_NUMBER: _ClassVar[int]
    MUT_NUC_FIELD_NUMBER: _ClassVar[int]
    CHROMOSOME_FIELD_NUMBER: _ClassVar[int]
    position: int
    ref_nuc: int
    par_nuc: int
    mut_nuc: _containers.RepeatedScalarFieldContainer[int]
    chromosome: str
    def __init__(self, position: _Optional[int] = ..., ref_nuc: _Optional[int] = ..., par_nuc: _Optional[int] = ..., mut_nuc: _Optional[_Iterable[int]] = ..., chromosome: _Optional[str] = ...) -> None: ...

class mutation_list(_message.Message):
    __slots__ = ("mutation",)
    MUTATION_FIELD_NUMBER: _ClassVar[int]
    mutation: _containers.RepeatedCompositeFieldContainer[mut]
    def __init__(self, mutation: _Optional[_Iterable[_Union[mut, _Mapping]]] = ...) -> None: ...

class condensed_node(_message.Message):
    __slots__ = ("node_name", "condensed_leaves")
    NODE_NAME_FIELD_NUMBER: _ClassVar[int]
    CONDENSED_LEAVES_FIELD_NUMBER: _ClassVar[int]
    node_name: str
    condensed_leaves: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, node_name: _Optional[str] = ..., condensed_leaves: _Optional[_Iterable[str]] = ...) -> None: ...

class node_metadata(_message.Message):
    __slots__ = ("clade_annotations",)
    CLADE_ANNOTATIONS_FIELD_NUMBER: _ClassVar[int]
    clade_annotations: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, clade_annotations: _Optional[_Iterable[str]] = ...) -> None: ...

class data(_message.Message):
    __slots__ = ("newick", "node_mutations", "condensed_nodes", "metadata")
    NEWICK_FIELD_NUMBER: _ClassVar[int]
    NODE_MUTATIONS_FIELD_NUMBER: _ClassVar[int]
    CONDENSED_NODES_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    newick: str
    node_mutations: _containers.RepeatedCompositeFieldContainer[mutation_list]
    condensed_nodes: _containers.RepeatedCompositeFieldContainer[condensed_node]
    metadata: _containers.RepeatedCompositeFieldContainer[node_metadata]
    def __init__(self, newick: _Optional[str] = ..., node_mutations: _Optional[_Iterable[_Union[mutation_list, _Mapping]]] = ..., condensed_nodes: _Optional[_Iterable[_Union[condensed_node, _Mapping]]] = ..., metadata: _Optional[_Iterable[_Union[node_metadata, _Mapping]]] = ...) -> None: ...
