from .sequences import router as sequences_router
from .databases import router as databases_router
from .files import router as files_router
from .alignment import router as alignment_router
from .crispr import router as crispr_router
from .variants import router as variants_router
from .primers import router as primers_router
from .pathways import router as pathways_router
from .analysis import router as analysis_router

__all__ = [
    "sequences_router", "databases_router", "files_router",
    "alignment_router", "crispr_router", "variants_router",
    "primers_router", "pathways_router", "analysis_router",
]
