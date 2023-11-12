import numpy as np
import pytest
from staver.MaxLFQ import MaxLFQOptimizer


# Happy path tests with various realistic test values
@pytest.mark.parametrize(
    "test_id, input_matrix, expected",
    [
        (
            "HP-1",
            np.array([[1, 2], [3, 4]]),
            {"estimate": np.array([1.5, 2.5]), "annotation": ""},
        ),
        (
            "HP-2",
            np.array([[np.nan, 2], [3, 4]]),
            {"estimate": np.array([3.5, 2.5]), "annotation": "1;1"},
        ),
        (
            "HP-3",
            np.array([[1, np.nan], [np.nan, 4]]),
            {"estimate": np.array([1, 4]), "annotation": "1;2"},
        ),
        # Add more happy path test cases
    ],
)
def test_maxLFQOptimizer_happy_path(test_id, input_matrix, expected):
    optimizer = MaxLFQOptimizer()

    # Act
    result = optimizer.maxLFQ_fast(input_matrix)

    # Assert
    assert np.allclose(
        result["estimate"], expected["estimate"], equal_nan=True
    ), f"Test {test_id} failed."
    assert result["annotation"] == expected["annotation"], f"Test {test_id} failed."


# Edge cases
@pytest.mark.parametrize(
    "test_id, input_matrix, expected",
    [
        (
            "EC-1",
            np.array([[np.nan, np.nan], [np.nan, np.nan]]),
            {"estimate": np.nan, "annotation": "All NaN values"},
        ),
        ("EC-2", np.array([[]]), {"estimate": np.nan, "annotation": "Empty array"}),
        ("EC-3", np.array([[1]]), {"estimate": np.array([1]), "annotation": [""]}),
        # Add more edge case test cases
    ],
)
def test_maxLFQOptimizer_edge_cases(test_id, input_matrix, expected):
    optimizer = MaxLFQOptimizer()

    # Act
    result = optimizer.maxLFQ_fast(input_matrix)

    # Assert
    assert (
        np.isnan(result["estimate"]).all() == np.isnan(expected["estimate"]).all()
    ), f"Test {test_id} failed."
    assert result["annotation"] == expected["annotation"], f"Test {test_id} failed."


# Error cases
@pytest.mark.parametrize(
    "test_id, input_matrix, error_type",
    [
        # Assuming the function should raise errors for invalid inputs
        # Uncomment and adjust the following lines if the function is expected to raise errors
        # ("ER-1", "invalid input", TypeError),
        # ("ER-2", None, TypeError),
        # Add more error case test cases
    ],
)
def test_maxLFQOptimizer_error_cases(test_id, input_matrix, error_type):
    optimizer = MaxLFQOptimizer()

    # Act / Assert
    with pytest.raises(error_type):
        optimizer.maxLFQ_fast(input_matrix)
