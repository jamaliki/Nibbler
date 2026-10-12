"""Lazy, cached public data for the structure-preparation regressions."""

import os
from pathlib import Path

import pytest

from tools.fetch_preparation_data import DEFAULT_DESTINATION, PreparationCorpus


@pytest.fixture(scope="session")
def preparation_data():
    destination = Path(os.environ.get("NIBBLER_TEST_DATA", DEFAULT_DESTINATION))
    return PreparationCorpus(
        destination, offline=os.environ.get("NIBBLER_TEST_DATA_OFFLINE") == "1"
    )
