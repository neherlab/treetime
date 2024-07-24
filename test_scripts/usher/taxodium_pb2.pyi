from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class MutationList(_message.Message):
    __slots__ = ("mutation",)
    MUTATION_FIELD_NUMBER: _ClassVar[int]
    mutation: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, mutation: _Optional[_Iterable[int]] = ...) -> None: ...

class MetadataSingleValuePerNode(_message.Message):
    __slots__ = ("metadata_name", "metadata_title", "mapping", "node_values")
    METADATA_NAME_FIELD_NUMBER: _ClassVar[int]
    METADATA_TITLE_FIELD_NUMBER: _ClassVar[int]
    MAPPING_FIELD_NUMBER: _ClassVar[int]
    NODE_VALUES_FIELD_NUMBER: _ClassVar[int]
    metadata_name: str
    metadata_title: str
    mapping: _containers.RepeatedScalarFieldContainer[str]
    node_values: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, metadata_name: _Optional[str] = ..., metadata_title: _Optional[str] = ..., mapping: _Optional[_Iterable[str]] = ..., node_values: _Optional[_Iterable[int]] = ...) -> None: ...

class AllNodeData(_message.Message):
    __slots__ = ("names", "x", "y", "countries", "dates", "lineages", "mutations", "parents", "genbanks", "epi_isl_numbers", "num_tips", "metadata_singles")
    NAMES_FIELD_NUMBER: _ClassVar[int]
    X_FIELD_NUMBER: _ClassVar[int]
    Y_FIELD_NUMBER: _ClassVar[int]
    COUNTRIES_FIELD_NUMBER: _ClassVar[int]
    DATES_FIELD_NUMBER: _ClassVar[int]
    LINEAGES_FIELD_NUMBER: _ClassVar[int]
    MUTATIONS_FIELD_NUMBER: _ClassVar[int]
    PARENTS_FIELD_NUMBER: _ClassVar[int]
    GENBANKS_FIELD_NUMBER: _ClassVar[int]
    EPI_ISL_NUMBERS_FIELD_NUMBER: _ClassVar[int]
    NUM_TIPS_FIELD_NUMBER: _ClassVar[int]
    METADATA_SINGLES_FIELD_NUMBER: _ClassVar[int]
    names: _containers.RepeatedScalarFieldContainer[str]
    x: _containers.RepeatedScalarFieldContainer[float]
    y: _containers.RepeatedScalarFieldContainer[float]
    countries: _containers.RepeatedScalarFieldContainer[int]
    dates: _containers.RepeatedScalarFieldContainer[int]
    lineages: _containers.RepeatedScalarFieldContainer[int]
    mutations: _containers.RepeatedCompositeFieldContainer[MutationList]
    parents: _containers.RepeatedScalarFieldContainer[int]
    genbanks: _containers.RepeatedScalarFieldContainer[str]
    epi_isl_numbers: _containers.RepeatedScalarFieldContainer[int]
    num_tips: _containers.RepeatedScalarFieldContainer[int]
    metadata_singles: _containers.RepeatedCompositeFieldContainer[MetadataSingleValuePerNode]
    def __init__(self, names: _Optional[_Iterable[str]] = ..., x: _Optional[_Iterable[float]] = ..., y: _Optional[_Iterable[float]] = ..., countries: _Optional[_Iterable[int]] = ..., dates: _Optional[_Iterable[int]] = ..., lineages: _Optional[_Iterable[int]] = ..., mutations: _Optional[_Iterable[_Union[MutationList, _Mapping]]] = ..., parents: _Optional[_Iterable[int]] = ..., genbanks: _Optional[_Iterable[str]] = ..., epi_isl_numbers: _Optional[_Iterable[int]] = ..., num_tips: _Optional[_Iterable[int]] = ..., metadata_singles: _Optional[_Iterable[_Union[MetadataSingleValuePerNode, _Mapping]]] = ...) -> None: ...

class AllData(_message.Message):
    __slots__ = ("node_data", "country_mapping", "lineage_mapping", "mutation_mapping", "date_mapping", "tree_description", "tree_title")
    NODE_DATA_FIELD_NUMBER: _ClassVar[int]
    COUNTRY_MAPPING_FIELD_NUMBER: _ClassVar[int]
    LINEAGE_MAPPING_FIELD_NUMBER: _ClassVar[int]
    MUTATION_MAPPING_FIELD_NUMBER: _ClassVar[int]
    DATE_MAPPING_FIELD_NUMBER: _ClassVar[int]
    TREE_DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    TREE_TITLE_FIELD_NUMBER: _ClassVar[int]
    node_data: AllNodeData
    country_mapping: _containers.RepeatedScalarFieldContainer[str]
    lineage_mapping: _containers.RepeatedScalarFieldContainer[str]
    mutation_mapping: _containers.RepeatedScalarFieldContainer[str]
    date_mapping: _containers.RepeatedScalarFieldContainer[str]
    tree_description: str
    tree_title: str
    def __init__(self, node_data: _Optional[_Union[AllNodeData, _Mapping]] = ..., country_mapping: _Optional[_Iterable[str]] = ..., lineage_mapping: _Optional[_Iterable[str]] = ..., mutation_mapping: _Optional[_Iterable[str]] = ..., date_mapping: _Optional[_Iterable[str]] = ..., tree_description: _Optional[str] = ..., tree_title: _Optional[str] = ...) -> None: ...
