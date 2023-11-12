import os
import pytest
import pandas as pd
from staver.data import read_proteomics_file

# Define file paths for testing
CSV_FILE = "/path/to/test.csv"
TSV_FILE = "/path/to/test.tsv"
XLS_FILE = "/path/to/test.xls"
MZID_FILE = "/path/to/test.mzid"
PEPXML_FILE = "/path/to/test.pepxml"
PROTXML_FILE = "/path/to/test.protxml"
MZTAB_FILE = "/path/to/test.mztab"
NON_EXISTENT_FILE = "/path/to/nonexistent.file"
UNSUPPORTED_FILE = "/path/to/unsupported.filetype"


# Create mock data for each file type
@pytest.fixture(scope="module", autouse=True)
def setup_files(tmpdir_factory):
    # Arrange
    tmpdir = tmpdir_factory.mktemp("data")
    csv_file = tmpdir.join("test.csv")
    tsv_file = tmpdir.join("test.tsv")
    xls_file = tmpdir.join("test.xls")
    mzid_file = tmpdir.join("test.mzid")
    pepxml_file = tmpdir.join("test.pepxml")
    protxml_file = tmpdir.join("test.protxml")
    mztab_file = tmpdir.join("test.mztab")
    unsupported_file = tmpdir.join("unsupported.filetype")

    # Create mock files with test data
    pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}).to_csv(csv_file, index=False)
    pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}).to_csv(
        tsv_file, sep="\t", index=False
    )
    pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}).to_excel(xls_file, index=False)
    # For specialized formats, create empty files as placeholders
    mzid_file.write("")
    pepxml_file.write("")
    protxml_file.write("")
    mztab_file.write("")
    unsupported_file.write("")

    # Return file paths
    return {
        "csv": str(csv_file),
        "tsv": str(tsv_file),
        "xls": str(xls_file),
        "mzid": str(mzid_file),
        "pepxml": str(pepxml_file),
        "protxml": str(protxml_file),
        "mztab": str(mztab_file),
        "unsupported": str(unsupported_file),
        "nonexistent": str(NON_EXISTENT_FILE),
    }


@pytest.mark.parametrize(
    "file_path, usecols, expected_columns, test_id",
    [
        # Happy path tests
        ("csv", ["col1", "col2"], ["col1", "col2"], "happy_csv"),
        ("tsv", ["col1"], ["col1"], "happy_tsv"),
        ("xls", ["col2"], ["col2"], "happy_xls"),
        # Edge cases
        ("csv", None, ["col1", "col2"], "edge_csv_no_usecols"),
        # Error cases
        ("nonexistent", ["col1", "col2"], ValueError, "error_nonexistent_file"),
        ("unsupported", ["col1", "col2"], ValueError, "error_unsupported_extension"),
    ],
)
def test_read_proteomics_file(
    setup_files, file_path, usecols, expected_columns, test_id
):
    file_path = setup_files[file_path]

    if test_id.startswith("happy") or test_id.startswith("edge"):
        # Act
        result = read_proteomics_file(file_path, usecols)

        # Assert
        assert list(result.columns) == expected_columns, f"Test failed for {test_id}"
    elif test_id.startswith("error"):
        # Act & Assert
        with pytest.raises(expected_columns) as exc_info:
            read_proteomics_file(file_path, usecols)
        assert str(file_path) in str(exc_info.value), f"Test failed for {test_id}"
