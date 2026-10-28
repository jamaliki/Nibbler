"""Skip tests that need Reduce3 chem_data, unless the run requires chem_data.

CI fetches the pinned subset (``python -m tools.fetch_chem_data``) and sets
``NIBBLER_REQUIRE_CHEM_DATA=1``, so a test that cannot find chem_data fails there
instead of skipping.
"""

from __future__ import annotations

import os
from typing import NoReturn

import pytest

REQUIRED = os.environ.get("NIBBLER_REQUIRE_CHEM_DATA") == "1"

needs_reduce3_chem_data = pytest.mark.skipif(
    not REQUIRED and not os.environ.get("REDUCE3_CHEM_DATA"),
    reason="Set REDUCE3_CHEM_DATA",
)


def skip_without_chem_data(reason: str) -> NoReturn:
    """Skip the test, or fail it when ``NIBBLER_REQUIRE_CHEM_DATA=1``."""
    if REQUIRED:
        pytest.fail(f"NIBBLER_REQUIRE_CHEM_DATA=1, but: {reason}", pytrace=False)
    pytest.skip(reason)
