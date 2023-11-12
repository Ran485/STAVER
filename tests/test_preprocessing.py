import pytest
import pandas as pd
import numpy as np
from staver.preprocessing import coefficient_variation, reduce_mem_usage


# Happy path tests with various realistic test values
@pytest.mark.parametrize(
    "data, log_transformed, dropna_thresh, expected",
    [
        # ID: Test with no NaN values and no log transformation
        (
            pd.DataFrame({"Protein1": [1, 2, 3], "Protein2": [4, 5, 6]}),
            False,
            0.5,
            pd.DataFrame(
                {
                    "Protein1": [1, 2, 3],
                    "Protein2": [4, 5, 6],
                    "Coefficient of Variation [%]": [50.0, 25.0, 16.666666666666664],
                    "Abundance": [2.5, 3.5, 4.5],
                    "Abundance_Rank": [3, 2, 1],
                }
            ),
        ),
        # ID: Test with NaN values and no log transformation
        (
            pd.DataFrame({"Protein1": [1, 2, np.nan], "Protein2": [4, np.nan, 6]}),
            False,
            0.5,
            pd.DataFrame(
                {
                    "Protein1": [1, 2, np.nan],
                    "Protein2": [4, np.nan, 6],
                    "Coefficient of Variation [%]": [100.0, np.nan, 0.0],
                    "Abundance": [2.5, np.nan, 6.0],
                    "Abundance_Rank": [2, np.nan, 1],
                }
            ).dropna(),
        ),
        # ID: Test with log transformation
        (
            pd.DataFrame({"Protein1": [1, 2, 4], "Protein2": [2, 4, 8]}),
            True,
            0.5,
            pd.DataFrame(
                {
                    "Protein1": [1, 2, 4],
                    "Protein2": [2, 4, 8],
                    "Coefficient of Variation [%]": [
                        33.33333333333333,
                        33.33333333333333,
                        33.33333333333333,
                    ],
                    "Abundance": [1.0, 2.0, 3.0],
                    "Abundance_Rank": [3, 2, 1],
                }
            ),
        ),
        # ID: Test with high dropna threshold
        (
            pd.DataFrame(
                {"Protein1": [1, np.nan, np.nan], "Protein2": [np.nan, np.nan, 2]}
            ),
            False,
            0.9,
            pd.DataFrame(
                {
                    "Protein1": [np.nan, np.nan, np.nan],
                    "Protein2": [np.nan, np.nan, 2],
                    "Coefficient of Variation [%]": [np.nan, np.nan, np.nan],
                    "Abundance": [np.nan, np.nan, 2.0],
                    "Abundance_Rank": [np.nan, np.nan, 1],
                }
            ).dropna(),
        ),
    ],
)
def test_happy_paths(data, log_transformed, dropna_thresh, expected):
    # Act
    result = coefficient_variation(data, log_transformed, dropna_thresh)

    # Assert
    pd.testing.assert_frame_equal(result, expected)


# Edge cases
@pytest.mark.parametrize(
    "data, log_transformed, dropna_thresh, expected",
    [
        # ID: Test with all NaN values
        (
            pd.DataFrame({"Protein1": [np.nan, np.nan], "Protein2": [np.nan, np.nan]}),
            False,
            0.5,
            pd.DataFrame(
                columns=[
                    "Protein1",
                    "Protein2",
                    "Coefficient of Variation [%]",
                    "Abundance",
                    "Abundance_Rank",
                ]
            ),
        ),
        # ID: Test with all zero values
        (
            pd.DataFrame({"Protein1": [0, 0], "Protein2": [0, 0]}),
            False,
            0.5,
            pd.DataFrame(
                {
                    "Protein1": [0, 0],
                    "Protein2": [0, 0],
                    "Coefficient of Variation [%]": [np.nan, np.nan],
                    "Abundance": [0.0, 0.0],
                    "Abundance_Rank": [1, 2],
                }
            ),
        ),
    ],
)
def test_edge_cases(data, log_transformed, dropna_thresh, expected):
    # Act
    result = coefficient_variation(data, log_transformed, dropna_thresh)

    # Assert
    pd.testing.assert_frame_equal(result, expected)


# Error cases
@pytest.mark.parametrize(
    "data, log_transformed, dropna_thresh, error",
    [
        # ID: Test with non-DataFrame input
        (None, False, 0.5, TypeError),
        # ID: Test with invalid dropna_thresh
        (pd.DataFrame({"Protein1": [1, 2], "Protein2": [3, 4]}), False, -1, ValueError),
        # ID: Test with invalid log_transformed type
        (
            pd.DataFrame({"Protein1": [1, 2], "Protein2": [3, 4]}),
            "not_bool",
            0.5,
            TypeError,
        ),
    ],
)
def test_error_cases(data, log_transformed, dropna_thresh, error):
    # Act and Assert
    with pytest.raises(error):
        coefficient_variation(data, log_transformed, dropna_thresh)


# test reduce_mem_usage
# Happy path tests with various realistic test values
@pytest.mark.parametrize(
    "data, expected_dtypes, expected_memory_decrease",
    [
        (
            {"int_column": [1, 2, 3], "float_column": [0.1, 0.2, 0.3]},
            {"int_column": "int8", "float_column": "float16"},
            True,
        ),
        (
            {"int_column": [1000, 2000, 3000], "float_column": [0.1, 0.2, 0.3]},
            {"int_column": "int16", "float_column": "float16"},
            True,
        ),
        (
            {"int_column": [1, 2, 3], "float_column": [1e-10, 2e-10, 3e-10]},
            {"int_column": "int8", "float_column": "float32"},
            True,
        ),
    ],
    ids=["small-int-float16", "medium-int-float16", "small-int-small-float"],
)
def test_reduce_mem_usage_happy_path(data, expected_dtypes, expected_memory_decrease):
    # Arrange
    df = pd.DataFrame(data)
    original_memory = df.memory_usage().sum()

    # Act
    df_reduced = reduce_mem_usage(df, verbose=False)
    reduced_memory = df_reduced.memory_usage().sum()

    # Assert
    for col, dtype in expected_dtypes.items():
        assert df_reduced[col].dtype == dtype
    assert (reduced_memory < original_memory) == expected_memory_decrease


# Edge cases
@pytest.mark.parametrize(
    "data, expected_dtypes",
    [
        ({"int_column": [0, 0, 0]}, {"int_column": "int8"}),
        ({"float_column": [np.nan, np.nan, np.nan]}, {"float_column": "float16"}),
    ],
    ids=["zero-int-column", "nan-float-column"],
)
def test_reduce_mem_usage_edge_cases(data, expected_dtypes):
    # Arrange
    df = pd.DataFrame(data)

    # Act
    df_reduced = reduce_mem_usage(df, verbose=False)

    # Assert
    for col, dtype in expected_dtypes.items():
        assert df_reduced[col].dtype == dtype


# Error cases
@pytest.mark.parametrize(
    "data, error_message",
    [
        (
            {"int_column": ["a", "b", "c"]},
            "Cannot convert non-numeric data to a numeric dtype.",
        ),
        (
            {"date_column": pd.date_range(start="1/1/2021", periods=3)},
            "Cannot downcast non-numeric dtype.",
        ),
    ],
    ids=["non-numeric-data", "non-numeric-dtype"],
)
def test_reduce_mem_usage_error_cases(data, error_message):
    # Arrange
    df = pd.DataFrame(data)

    # Act / Assert
    with pytest.raises(TypeError) as exc_info:
        reduce_mem_usage(df, verbose=False)
    assert error_message in str(exc_info.value)
