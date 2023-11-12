import pytest
from unittest.mock import patch, MagicMock
from staver.data import joblib_load_file
from pandas.testing import assert_frame_equal
import pandas as pd

# Assuming NUM_WORKERS is defined somewhere in the module
from staver.data import NUM_WORKERS


# Mocking external dependencies
def mock_get_all_files(inpath, extension):
    return [f"{inpath}/file1{extension}", f"{inpath}/file2{extension}"]


def mock_read_file(file):
    return pd.DataFrame([[1, 2], [3, 4]], columns=["Column1", "Column2"])


@pytest.mark.parametrize(
    "test_id, inpath, extension, num_workers, expected",
    [
        # Happy path tests
        (
            "HP-1",
            "/path/to/files",
            ".tsv",
            4,
            pd.DataFrame(
                [[1, 2], [3, 4], [1, 2], [3, 4]], columns=["Column1", "Column2"]
            ),
        ),
        (
            "HP-2",
            "/path/to/files",
            ".csv",
            2,
            pd.DataFrame(
                [[1, 2], [3, 4], [1, 2], [3, 4]], columns=["Column1", "Column2"]
            ),
        ),
        # Edge cases
        (
            "EC-1",
            "/path/to/files",
            ".tsv",
            1,
            pd.DataFrame([[1, 2], [3, 4]], columns=["Column1", "Column2"]),
        ),  # Single worker
        (
            "EC-2",
            "/path/to/files",
            ".tsv",
            0,
            pd.DataFrame([], columns=["Column1", "Column2"]),
        ),  # Zero workers
        # Error cases
        ("ER-1", None, ".tsv", 4, ValueError),  # inpath is None
        ("ER-2", "/path/to/files", None, 4, ValueError),  # extension is None
        (
            "ER-3",
            "/path/to/files",
            ".tsv",
            -1,
            ValueError,
        ),  # Negative number of workers
    ],
)
def test_joblib_load_file(test_id, inpath, extension, num_workers, expected):
    # Arrange
    with patch("staver.data.get_all_files", side_effect=mock_get_all_files), patch(
        "staver.data.read_file", side_effect=mock_read_file
    ), patch("staver.data.ProcessPoolExecutor") as mock_executor, patch(
        "staver.data.Progress"
    ) as mock_progress:
        mock_executor.return_value.__enter__.return_value = mock_executor
        mock_executor.submit.side_effect = lambda f, file: MagicMock(
            result=lambda: mock_read_file(file)
        )
        mock_progress.return_value.__enter__.return_value = mock_progress
        mock_progress.add_task.return_value = 1

        # Act
        if isinstance(expected, pd.DataFrame):
            result = joblib_load_file(inpath, extension, num_workers)

        # Assert
        if isinstance(expected, pd.DataFrame):
            assert_frame_equal(result, expected)
        else:
            with pytest.raises(expected):
                joblib_load_file(inpath, extension, num_workers)
