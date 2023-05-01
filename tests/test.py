import pandas as pd
import sys

sys.path.append("staver")

from utils import reduce_mem_usage


def test_reduce_mem_usage():
    # Create a DataFrame with columns of different data types
    df = pd.DataFrame(
        {
            'int_col': [1, 2, 3],
            'float_col': [1.0, 2.0, 3.0],
            'str_col': ['foo', 'bar', 'baz'],
        }
    )

    # Call the function to reduce memory usage
    reduced_df = reduce_mem_usage(df)

    # Check that the function returns a DataFrame
    assert isinstance(reduced_df, pd.DataFrame)

    # Check that the numeric columns have been converted to smaller data types
    assert str(reduced_df['int_col'].dtype) == 'int8'
    assert str(reduced_df['float_col'].dtype) == 'float16'

    # Check that the string column is unchanged
    assert str(reduced_df['str_col'].dtype) == 'object'

    # Check that the function prints memory usage information if verbose is True
    assert reduce_mem_usage(df, verbose=True) == reduced_df
