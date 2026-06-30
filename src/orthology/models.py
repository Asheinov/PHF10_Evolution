"""
Центральные структуры данных ортологического анализа.

Содержит:

- OrthologHit
    Одна запись из конкретного источника
    (OMA, OrthoDB, Ensembl, OrthoFinder).

- OrthologCandidate
    Объединённый кандидат после агрегации
    нескольких источников.

- OrthologSet
    Финальный набор ортологов,
    разделённый по уровням уверенности.

Этот модуль не содержит бизнес-логики.
Только описание объектов данных.
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict


# --------------------------------------------------
# Исходная запись из одной базы ортологии
# --------------------------------------------------

@dataclass(slots=True)
class OrthologHit:
    """
    Одна запись, полученная из конкретного источника.

    Пример:
        PHF10 (Homo sapiens)
            ->
        PHF10 (Mus musculus)

    из OMA или Ensembl.
    """

    source: str

    query_gene: str
    target_gene: str

    query_species: str
    target_species: str

    query_protein_id: Optional[str] = None
    target_protein_id: Optional[str] = None

    orthology_type: Optional[str] = None
    # Например:
    # 1:1
    # 1:many
    # many:many

    raw_score: Optional[float] = None

    tree_support: Optional[float] = None

    metadata: Dict[str, str] = field(default_factory=dict)


# --------------------------------------------------
# Объединённый кандидат после агрегации
# --------------------------------------------------

@dataclass(slots=True)
class OrthologCandidate:
    """
    Кандидат в ортологи после объединения
    данных из нескольких источников.

    Для одного target_gene может существовать
    несколько подтверждений из разных баз.
    """

    target_gene: str
    target_species: str

    hits: List[OrthologHit] = field(default_factory=list)

    sources: List[str] = field(default_factory=list)

    n_sources: int = 0

    mean_score: Optional[float] = None
    max_score: Optional[float] = None
    combined_score: Optional[float] = None

    confidence_class: Optional[str] = None
    # HIGH
    # MEDIUM
    # LOW

    conflict_flag: bool = False

    metadata: Dict[str, str] = field(default_factory=dict)


# --------------------------------------------------
# Финальный набор ортологов
# --------------------------------------------------

@dataclass(slots=True)
class OrthologSet:
    """
    Финальный набор ортологов после
    агрегации и классификации.

    Используется как основной объект,
    который передаётся в последующие
    этапы анализа.
    """

    query_gene: str

    high_confidence: List[OrthologCandidate] = field(
        default_factory=list
    )

    medium_confidence: List[OrthologCandidate] = field(
        default_factory=list
    )

    low_confidence: List[OrthologCandidate] = field(
        default_factory=list
    )

    excluded: List[OrthologCandidate] = field(
        default_factory=list
    )


# --------------------------------------------------
# Константы уровней уверенности
# --------------------------------------------------

CONFIDENCE_HIGH = "HIGH"
CONFIDENCE_MEDIUM = "MEDIUM"
CONFIDENCE_LOW = "LOW"


# --------------------------------------------------
# Константы типов ортологии
# --------------------------------------------------

ORTHOLOGY_ONE_TO_ONE = "1:1"
ORTHOLOGY_ONE_TO_MANY = "1:many"
ORTHOLOGY_MANY_TO_MANY = "many:many"
