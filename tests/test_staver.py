import os
import pytest
from unittest.mock import patch, MagicMock
from staver import Staver, StaverConfig

# Constants for tests
VALID_INPATH = "/valid/inpath/"
VALID_OUTPATH = "/valid/outpath/"
VALID_SUFFIX = "_test_suffix"
VALID_THRESHOLD_NUMS = 4
VALID_WITH_STANDARD_DATA = True
VALID_DIA_PATH = "/valid/dia_path/"
VALID_STANDARD_DATASET_PATH = "/valid/standard_dataset_path/"
VALID_DIA_PEP_DATA_OUTPATH = "/valid/dia_pep_data_outpath/"
VALID_DIA_PROTEIN_DATA_OUTPATH = "/valid/dia_protein_data_outpath/"
VALID_FDR_THRESHOLD = 0.01
VALID_COUNT_CUTOFF_SAME_FILES = 1
VALID_COUNT_CUTOFF_DIFF_FILES = 2
VALID_PEPTIDES_CV_THRESH = 0.3
VALID_PROTEINS_CV_THRESH = 0.3
VALID_TOP_PEPTIDES_NUMS = 6
VALID_NA_THRESHOLD = 0.3


# Helper function to create a valid StaverConfig instance
def create_valid_config():
    return StaverConfig(
        threshold_nums=VALID_THRESHOLD_NUMS,
        with_standard_data=VALID_WITH_STANDARD_DATA,
        dia_path=VALID_DIA_PATH,
        standard_dataset_path=VALID_STANDARD_DATASET_PATH,
        dia_peptide_data_outpath=VALID_DIA_PEP_DATA_OUTPATH,
        dia_protein_data_outpath=VALID_DIA_PROTEIN_DATA_OUTPATH,
        fdr_threshold=VALID_FDR_THRESHOLD,
        count_cutoff_same_files=VALID_COUNT_CUTOFF_SAME_FILES,
        count_cutoff_diff_files=VALID_COUNT_CUTOFF_DIFF_FILES,
        peptides_cv_thresh=VALID_PEPTIDES_CV_THRESH,
        proteins_cv_thresh=VALID_PROTEINS_CV_THRESH,
        top_peptides_nums=VALID_TOP_PEPTIDES_NUMS,
        na_threshold=VALID_NA_THRESHOLD,
        file_suffix=VALID_SUFFIX,
    )


# Parametrized test for StaverConfig initialization
@pytest.mark.parametrize(
    "threshold_nums, expected",
    [
        (
            None,
            os.cpu_count() - 2 if os.cpu_count() else 2,
        ),  # ID: CONFIG-THRESHOLD-DEFAULT
        (4, 4),  # ID: CONFIG-THRESHOLD-VALID
        (-1, ValueError),  # ID: CONFIG-THRESHOLD-INVALID
    ],
)
def test_staver_config_initialization(threshold_nums, expected):
    # Arrange
    if expected is ValueError:
        with pytest.raises(ValueError):
            StaverConfig(threshold_nums=threshold_nums)
    else:
        # Act
        config = StaverConfig(threshold_nums=threshold_nums)
        # Assert
        assert config.threshold_nums == expected


# Parametrized test for Staver.preprocess_data
@pytest.mark.parametrize(
    "inpath, file_suffix, expected_exception",
    [
        (VALID_INPATH, VALID_SUFFIX, None),  # ID: PREPROCESS-HAPPY-PATH
        ("/invalid/path/", VALID_SUFFIX, Exception),  # ID: PREPROCESS-INVALID-INPATH
        (VALID_INPATH, None, Exception),  # ID: PREPROCESS-INVALID-SUFFIX
    ],
)
def test_preprocess_data(inpath, file_suffix, expected_exception):
    # Arrange
    staver = Staver()
    if expected_exception:
        with pytest.raises(expected_exception):
            # Act
            staver.preprocess_data(inpath, file_suffix)
    else:
        with patch("os.makedirs") as mock_makedirs, patch(
            "os.path.exists", return_value=False
        ), patch("staver.merge_data") as mock_merge_data:
            # Act
            staver.preprocess_data(inpath, file_suffix)
            # Assert
            mock_makedirs.assert_called_once()
            mock_merge_data.assert_called_once()


# Parametrized test for Staver.run_staver
@pytest.mark.parametrize(
    "with_standard_data, expected_exception",
    [
        (True, None),  # ID: RUN-STANDARD-HAPPY-PATH
        (False, None),  # ID: RUN-NOSTANDARD-HAPPY-PATH
        (True, Exception),  # ID: RUN-STANDARD-EXCEPTION
        (False, Exception),  # ID: RUN-NOSTANDARD-EXCEPTION
    ],
)
def test_run_staver(with_standard_data, expected_exception):
    # Arrange
    config = create_valid_config()
    config.with_standard_data = with_standard_data
    staver = Staver()
    if expected_exception:
        with pytest.raises(expected_exception), patch.object(
            staver.highCIPeptides, "main", side_effect=Exception
        ), patch.object(
            staver.peptideToProteinInference, "parallel_main", side_effect=Exception
        ):
            # Act
            staver.run_staver(
                number_threshods=config.threshold_nums,
                with_standard_data=config.with_standard_data,
                dia_path=config.dia_path,
                dia_pep_data_outpath=config.dia_peptide_data_outpath,
                standard_dataset_path=config.standard_dataset_path,
                file_suffix=config.file_suffix,
                count_cutoff_same_files=config.count_cutoff_same_files,
                count_cutoff_diff_files=config.count_cutoff_diff_files,
                peptides_cv_thresh=config.peptides_cv_thresh,
                dia_protein_data_outpath=config.dia_protein_data_outpath,
                proteins_cv_thresh=config.proteins_cv_thresh,
                top_peptides_nums=config.top_peptides_nums,
                fdr_threshold=config.fdr_threshold,
                na_threshold=config.na_threshold,
            )
    else:
        with patch.object(
            staver.highCIPeptides, "main"
        ) as mock_highCIPeptides_main, patch.object(
            staver.peptideToProteinInference, "parallel_main"
        ) as mock_peptideToProteinInference_parallel_main, patch.object(
            staver, "preprocess_data"
        ) as mock_preprocess_data:
            # Act
            staver.run_staver(
                number_threshods=config.threshold_nums,
                with_standard_data=config.with_standard_data,
                dia_path=config.dia_path,
                dia_pep_data_outpath=config.dia_peptide_data_outpath,
                standard_dataset_path=config.standard_dataset_path,
                file_suffix=config.file_suffix,
                count_cutoff_same_files=config.count_cutoff_same_files,
                count_cutoff_diff_files=config.count_cutoff_diff_files,
                peptides_cv_thresh=config.peptides_cv_thresh,
                dia_protein_data_outpath=config.dia_protein_data_outpath,
                proteins_cv_thresh=config.proteins_cv_thresh,
                top_peptides_nums=config.top_peptides_nums,
                fdr_threshold=config.fdr_threshold,
                na_threshold=config.na_threshold,
            )
            # Assert
            mock_highCIPeptides_main.assert_called_once()
            mock_peptideToProteinInference_parallel_main.assert_called_once()
            mock_preprocess_data.assert_called_once()
