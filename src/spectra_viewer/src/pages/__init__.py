from .page1 import Page1
from .peak_list import peak_list
from .mgf_page import mgf_page
from .compare import compare
from ..utils import Page

from typing import Dict, Type


PAGE_MAP: Dict[str, Type[Page]] = {
    "MGF File": mgf_page,
    "Peak List": peak_list,
    # "UNPD-ISDB": Page1,
}

PAGE_COMP: Dict[str, Type[Page]] = {
    "Compare": compare
}

__all__ = ["PAGE_MAP"]