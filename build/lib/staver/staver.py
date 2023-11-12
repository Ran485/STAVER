import os
from HighCIPeptides import HighCIPeptides
from PeptideToProteinInference import PeptideToProteinInference
from merge_data import merge_data
from utils import create_output_dir, timer


class Staver:
    """
    Class for running the STAVER algorithm for proteomics data analysis.

    This class encapsulates the workflow of the STAVER algorithm, which involves
    preprocessing data, peptide to protein inference, and other steps necessary
    for proteomics data analysis.

    Example:
        config = StaverConfig(
            threshold_nums=16,
            with_standard_data=True,
            dia_path="~/data/dia_raw/",
            standard_dataset_path="~/data/dia_standard_dataset/",
            # ... other configurations ...
        )
        staver = Staver()
        staver.run_staver(
            number_threshods=config.threshold_nums,
            with_standard_data=config.with_standard_data,
            dia_path=config.dia_path,
            dia_pep_data_outpath=config.dia_peptide_data_outpath,
            # ... other parameters ...
        )

    Attributes:
        highCIPeptides (HighCIPeptides): Instance for high confidence peptides processing.
        peptideToProteinInference (PeptideToProteinInference): Instance for peptide to protein inference.
    """

    def __init__(self):
        """
        Initializes an instance of the Staver class.
        """
        self.highCIPeptides = HighCIPeptides()
        self.peptideToProteinInference = PeptideToProteinInference()

    def preprocess_data(self, inpath: str, file_suffix) -> None:
        """Preprocess the data.

        Preprocesses the data by creating a processed data directory, merging the input data,
        and storing the preprocessed data in the directory.

        Args:
            inpath (str): The input data path.
            file_suffix: The suffix for the processed file.

        Returns:
            None

        Raises:
            Exception: If an error occurs during the preprocessing of the data.
        """
        try:
            outpath = os.path.join(
                os.path.dirname(os.path.dirname(inpath)), "processed_data/"
            )
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            print(f"\nPreprocessed data stored in `{outpath}`\n")
            file_extend = f"{file_suffix}.csv"
            merge_data(
                inpath,
                usecols=["Gene_Symbol", "median"],
                filetype_extension=file_extend,
            ).to_csv(f"{outpath}Protein_matrix.csv", index=False)

        except Exception as e:
            print(f"\nWhile merge`{inpath}`\n, There occurred an error:\n + {str(e)}\n")

    @timer
    def run_staver(
        self,
        number_threshods: int,
        with_standard_data: bool,
        dia_path: str,
        standard_dataset_path: str,
        dia_pep_data_outpath: str,
        dia_protein_data_outpath: str,
        file_suffix: str,
        fdr_threshold: float,
        count_cutoff_same_files: int,
        count_cutoff_diff_files: int,
        peptides_cv_thresh: float,
        proteins_cv_thresh: float,
        top_peptides_nums: int,
        na_threshold: float,
    ) -> None:
        """
        Main function to run the STAVER algorithm.

        Executes the various steps of the STAVER algorithm, including data preprocessing,
        peptide to protein inference, and additional analysis steps.

        Args:
            number_threshods (int): Number of thresholds for computational operations.
            with_standard_data (bool): Flag to indicate usage of standard reference data.
            dia_path (str): Path to the DIA data directory.
            standard_dataset_path (str): Path to the standard reference dataset.
            dia_pep_data_outpath (str): Output path for processed peptide data.
            dia_protein_data_outpath (str): Output path for processed protein data.
            file_suffix (str): Suffix for file identification.
            fdr_threshold (float): Threshold for FDR.
            count_cutoff_same_files (int): Cutoff count for same files.
            count_cutoff_diff_files (int): Cutoff count for different files.
            peptides_cv_thresh (float): CV threshold for peptides.
            proteins_cv_thresh (float): CV threshold for proteins.
            top_peptides_nums (int): Number of top precursor ions.
            na_threshold (float): Threshold for null peptides.

        Returns:
            None

        Raises:
            Exception: If an error occurs during the execution of the STAVER algorithm.
        """
        if with_standard_data:
            try:
                self.highCIPeptides.num_workers = number_threshods
                self.highCIPeptides.reference_dataset_path = standard_dataset_path
                self.highCIPeptides.count_cutoff_same_libs = count_cutoff_same_files
                self.highCIPeptides.count_cutoff_diff_libs = count_cutoff_diff_files
                self.highCIPeptides.cv_thresh_of_same_libs = peptides_cv_thresh
                self.highCIPeptides.cv_thresh_of_diff_libs = peptides_cv_thresh
                self.highCIPeptides.top_precursor_ions = top_peptides_nums
                self.highCIPeptides.fdr_threshold = fdr_threshold
                self.highCIPeptides.outpath = dia_pep_data_outpath
                reference_dataset = self.highCIPeptides.main()

                self.peptideToProteinInference.file_suffix = file_suffix
                self.peptideToProteinInference.input_data_path = dia_path
                self.peptideToProteinInference.reference_dataset = reference_dataset
                self.peptideToProteinInference.na_threshold = na_threshold
                self.peptideToProteinInference.cv_threshold = proteins_cv_thresh
                self.peptideToProteinInference.outpath = dia_protein_data_outpath
                self.peptideToProteinInference.parallel_main()

                self.preprocess_data(dia_protein_data_outpath, file_suffix)

            except Exception as e:
                error_message = (
                    "\nWhile processing `"
                    + dia_path
                    + " WITH_STANDARDE_DATA`\n, There occurred an error:\n"
                    + str(e)
                )
                print(error_message)
        else:
            try:
                self.highCIPeptides.num_workers = number_threshods
                print("Load the default standard dataset for calibration!")

                self.peptideToProteinInference.file_suffix = file_suffix
                self.peptideToProteinInference.input_data_path = dia_path
                self.peptideToProteinInference.outpath = dia_protein_data_outpath
                self.peptideToProteinInference.na_threshold = na_threshold
                self.peptideToProteinInference.cv_threshold = proteins_cv_thresh
                self.peptideToProteinInference.parallel_main()

                self.preprocess_data(dia_protein_data_outpath, file_suffix)

            except Exception as e:
                error_message = (
                    "\nWhile processing `"
                    + dia_path
                    + " WITHOUT_STANDARDE_DATA`\n, There occurred an error:\n"
                    + str(e)
                )
                print(error_message)


class StaverConfig:
    """
    Configuration class for the STAVER algorithm.

    Holds all the configuration settings required to run the STAVER algorithm.

    Attributes:
        threshold_nums (int): Number of thresholds for computational operations.
        with_standard_data (bool): Flag to indicate usage of standard reference data.
        # ... other attributes ...

    Example:
        config = StaverConfig(
            threshold_nums=16,
            with_standard_data=True,
            dia_path="~/data/dia_raw/",
            standard_dataset_path="~/data/dia_standard_dataset/",
            # ... other configurations ...
        )
    """

    def __init__(
        self,
        threshold_nums: int = None,
        with_standard_data: bool = True,
        count_cutoff_same_files: int = 1,
        count_cutoff_diff_files: int = 2,
        peptides_cv_thresh: float = 0.3,
        top_peptides_nums: int = 6,
        fdr_threshold: float = 0.01,
        dia_path: str = None,
        standard_dataset_path: str = "~/data/dian-raw/",
        dia_peptide_data_outpath: str = create_output_dir("results/peptides"),
        dia_protein_data_outpath: str = create_output_dir("results/proteins"),
        proteins_cv_thresh: float = 0.3,
        na_threshold: float = 0.3,
        file_suffix: str = "_F1_R1",
        **kwargs,
    ):
        """
        Initializes a new instance of the StaverConfig class.

        Args:
            threshold_nums (int): Number of thresholds for computational operations.
            with_standard_data (bool): Flag to indicate usage of standard reference data.
            # ... other parameters ...

        Returns:
            None
        """
        self.threshold_nums = self._validate_threshold_nums(threshold_nums)
        self.with_standard_data = with_standard_data
        self.count_cutoff_same_files = count_cutoff_same_files
        self.count_cutoff_diff_files = count_cutoff_diff_files
        self.peptides_cv_thresh = peptides_cv_thresh
        self.top_peptides_nums = top_peptides_nums
        self.fdr_threshold = fdr_threshold
        self.dia_path = dia_path
        self.standard_dataset_path = standard_dataset_path
        self.dia_peptide_data_outpath = dia_peptide_data_outpath
        self.dia_protein_data_outpath = dia_protein_data_outpath
        self.proteins_cv_thresh = proteins_cv_thresh
        self.na_threshold = na_threshold
        self.file_suffix = file_suffix
        # ... other parameters ...
        # ... initialize other attributes ...
        # Additional parameters passed through kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

    def _validate_threshold_nums(self, threshold_nums):
        """Validate or set default value for threshold_nums."""
        if threshold_nums is not None:
            if not isinstance(threshold_nums, int) or threshold_nums < 1:
                raise ValueError("threshold_nums must be a positive integer.")
            return threshold_nums
        else:
            cpu_count = os.cpu_count() or 4  # Default to 4 if cpu_count is None
            return max(1, cpu_count - 2)


if __name__ == "__main__":
    # Set the parameters for the StaverConfig instance config

    # The StaverConfig instance config is used to pass the configuration values to the `run_staver`
    # function. Each configuration value is accessed through the corresponding attribute of the config
    # object. This ensures that the function receives the correct configuration values for execution.
    config = StaverConfig()

    # Create an instance of the Staver class
    staver = Staver()

    # Run the staver algorithm
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
