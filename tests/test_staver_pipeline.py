import pytest
import os
from unittest.mock import patch, mock_open
from staver.staver_pipeline import (
    ValidateParser,
    print_arguments,
    export_arguments,
    dict_to_arglist,
    staver_pipeline,
)

# Constants for tests
VALID_INPUT_PATH = "/valid/input/path"
VALID_OUTPUT_PATH = "/valid/output/path"
INVALID_PATH = "/invalid/path"
NON_DIR_PATH = "/not/a/directory"
VALID_INT = 4
INVALID_INT = "four"
VALID_FLOAT = 0.5
INVALID_FLOAT = "half"
VALID_ARGS_DICT = {
    "thread_numbers": VALID_INT,
    "input": VALID_INPUT_PATH,
    "output_peptide": VALID_OUTPUT_PATH,
    "output_protein": VALID_OUTPUT_PATH,
    "fdr": VALID_FLOAT,
    "count_cutoff_same_libs": VALID_INT,
    "count_cutoff_diff_libs": VALID_INT,
    "peptides_cv_thresh": VALID_FLOAT,
    "proteins_cv_thresh": VALID_FLOAT,
    "na_threshold": VALID_FLOAT,
    "top_precursor_ions": VALID_INT,
    "normalization_method": "median",
    "file_suffix": "_F1_R1",
    "sample_type": "test_sample",
    "verbose": False,
}
VALID_ARGS_LIST = dict_to_arglist(VALID_ARGS_DICT)


# Helper function to create a temporary directory structure
@pytest.fixture
def create_temp_directory(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    return input_dir, output_dir


# Parametrized test for ValidateParser.valid_input_path
@pytest.mark.parametrize(
    "test_id, input_path, expected",
    [
        ("happy_path", VALID_INPUT_PATH, VALID_INPUT_PATH),
        ("error_nonexistent", INVALID_PATH, pytest.raises(argparse.ArgumentTypeError)),
        (
            "error_not_a_directory",
            NON_DIR_PATH,
            pytest.raises(argparse.ArgumentTypeError),
        ),
    ],
)
def test_valid_input_path(test_id, input_path, expected, create_temp_directory):
    # Arrange
    input_dir, _ = create_temp_directory
    if test_id == "happy_path":
        input_path = str(input_dir)

    # Act and Assert
    if isinstance(expected, type) and issubclass(expected, Exception):
        with expected:
            ValidateParser.valid_input_path(input_path)
    else:
        assert ValidateParser.valid_input_path(input_path) == expected


# Parametrized test for ValidateParser.valid_output_path
@pytest.mark.parametrize(
    "test_id, output_path, expected",
    [
        ("happy_path", VALID_OUTPUT_PATH, VALID_OUTPUT_PATH),
        (
            "error_not_a_directory",
            NON_DIR_PATH,
            pytest.raises(argparse.ArgumentTypeError),
        ),
    ],
)
def test_valid_output_path(test_id, output_path, expected, create_temp_directory):
    # Arrange
    _, output_dir = create_temp_directory
    if test_id == "happy_path":
        output_path = str(output_dir)

    # Act and Assert
    if isinstance(expected, type) and issubclass(expected, Exception):
        with expected:
            ValidateParser.valid_output_path(output_path)
    else:
        assert ValidateParser.valid_output_path(output_path) == expected


# Parametrized test for ValidateParser.valid_int
@pytest.mark.parametrize(
    "test_id, value, expected",
    [
        ("happy_path", str(VALID_INT), VALID_INT),
        ("error_non_int", INVALID_INT, pytest.raises(argparse.ArgumentTypeError)),
        ("error_negative", "-1", pytest.raises(argparse.ArgumentTypeError)),
    ],
)
def test_valid_int(test_id, value, expected):
    # Act and Assert
    if isinstance(expected, type) and issubclass(expected, Exception):
        with expected:
            ValidateParser.valid_int(value)
    else:
        assert ValidateParser.valid_int(value) == expected


# Parametrized test for ValidateParser.valid_float
@pytest.mark.parametrize(
    "test_id, value, expected",
    [
        ("happy_path", str(VALID_FLOAT), VALID_FLOAT),
        ("error_non_float", INVALID_FLOAT, pytest.raises(argparse.ArgumentTypeError)),
        ("error_out_of_range", "1.5", pytest.raises(argparse.ArgumentTypeError)),
    ],
)
def test_valid_float(test_id, value, expected):
    # Act and Assert
    if isinstance(expected, type) and issubclass(expected, Exception):
        with expected:
            ValidateParser.valid_float(value)
    else:
        assert ValidateParser.valid_float(value) == expected


# Parametrized test for print_arguments
@pytest.mark.parametrize(
    "test_id, args", [("happy_path", argparse.Namespace(**VALID_ARGS_DICT))]
)
def test_print_arguments(test_id, args, capsys):
    # Act
    print_arguments(args)

    # Assert
    captured = capsys.readouterr()
    for key, value in VALID_ARGS_DICT.items():
        assert f"{key}: {value}" in captured.out


# Parametrized test for export_arguments
@pytest.mark.parametrize(
    "test_id, args, output_path",
    [("happy_path", argparse.Namespace(**VALID_ARGS_DICT), VALID_OUTPUT_PATH)],
)
def test_export_arguments(test_id, args, output_path, create_temp_directory):
    # Arrange
    _, output_dir = create_temp_directory
    output_path = str(output_dir)
    expected_file_content = (
        "All parsed arguments:  \r "
        + "  \r ".join([f"{key}: {value}" for key, value in VALID_ARGS_DICT.items()])
        + "  \r "
    )

    # Act
    export_arguments(args, output_path)

    # Assert
    with open(os.path.join(output_path, "staver_arguments.txt"), "r") as f:
        file_content = f.read()
        assert file_content == expected_file_content


# Parametrized test for dict_to_arglist
@pytest.mark.parametrize(
    "test_id, arg_dict, expected", [("happy_path", VALID_ARGS_DICT, VALID_ARGS_LIST)]
)
def test_dict_to_arglist(test_id, arg_dict, expected):
    # Act
    result = dict_to_arglist(arg_dict)

    # Assert
    assert result == expected


# Parametrized test for staver_pipeline
@pytest.mark.parametrize(
    "test_id, args, expected_output, expected_exception",
    [
        ("happy_path", VALID_ARGS_LIST, "All parsed arguments:", None),
        ("error_invalid_input", ["--input", INVALID_PATH], None, SystemExit),
    ],
)
def test_staver_pipeline(
    test_id, args, expected_output, expected_exception, create_temp_directory, capsys
):
    # Arrange
    input_dir, output_dir = create_temp_directory
    args = [
        arg.replace(VALID_INPUT_PATH, str(input_dir)).replace(
            VALID_OUTPUT_PATH, str(output_dir)
        )
        for arg in args
    ]

    # Act
    if expected_exception:
        with pytest.raises(expected_exception):
            staver_pipeline(args)
    else:
        staver_pipeline(args)

        # Assert
        captured = capsys.readouterr()
        assert expected_output in captured.out
