import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from sklearn.cluster import DBSCAN
from pyod.models.knn import KNN
from pyod.models.iforest import IForest
from sklearn.covariance import EllipticEnvelope
from sklearn.model_selection import GridSearchCV
from joblib import Parallel, delayed

# python matplotlib export editable PDF
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["figure.dpi"] = 150


class OutlierDetection:
    """A class for outlier detection using various methods and visualization using Plotly.

    Attributes:
        contamination (float): The proportion of outliers in the data.
    """

    def __init__(self, contamination=0.05):
        """Initializes the OutlierDetectionPlotly with a contamination rate."""
        self.contamination = contamination

    def _ocsvm(self, X_scaled):
        """
        Helper method to fit and predict using OneClassSVM.

        Args:
            X_scaled (numpy.ndarray): The standardized data for outlier detection.

        Returns:
            numpy.ndarray: Array of outlier predictions (1: outlier, 0: inlier).
        """
        model = OneClassSVM(nu=self.contamination)
        model.fit(X_scaled)
        return model.predict(X_scaled)

    def optimize_ocsvm_parameters(self, X, param_grid=None):
        """
        Simplified method to optimize parameters for OneClassSVM.

        Args:
            X (numpy.ndarray): The data for outlier detection.
            param_grid (dict, optional): The parameter grid.

        Returns:
            OneClassSVM: The OneClassSVM model with selected parameters.
        """
        X_scaled = self.scaler.fit_transform(X)

        if param_grid is None:
            param_grid = {
                "nu": [0.01, 0.05, 0.1],
                "gamma": ["scale", "auto", 0.1, 0.01],
            }

        best_score = -np.inf
        best_model = None

        for nu in param_grid["nu"]:
            for gamma in param_grid["gamma"]:
                model = OneClassSVM(nu=nu, gamma=gamma)
                model.fit(X_scaled)
                score = model.score_samples(X_scaled).mean()  # Example scoring logic
                if score > best_score:
                    best_score = score
                    best_model = model

        return best_model

    def detect_outliers(self, X, method="iforest"):
        """Detects outliers in the given data using a specified method.

        Args:
            X (numpy.ndarray): The data for outlier detection.
            method (str): The method for outlier detection.

        Returns:
            numpy.ndarray: Predictions where 1 indicates an outlier and 0 an inlier.
        """
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Initialize the specified outlier detection model
        if method == "iforest":
            model = IForest(contamination=self.contamination, random_state=42)
        elif method == "knn":
            model = KNN(contamination=self.contamination)
        elif method == "lof":
            model = LocalOutlierFactor(n_neighbors=20, contamination=self.contamination)
        elif method == "ocsvm":
            model = self.optimize_ocsvm_parameters(X)
        elif method == "dbscan":
            model = DBSCAN(eps=3.0, min_samples=100)
        # elif method == 'autoencoder':
        #     model = AutoEncoder(hidden_neurons=[64, 32, 32, 64])
        elif method == "ee":
            model = EllipticEnvelope(contamination=self.contamination)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Fit model and predict outliers
        if method in ["iforest", "knn", "ocsvm", "ee"]:
            model.fit(X_scaled)
            predictions = model.predict(X_scaled)
            if method == "ocsvm":
                predictions = np.where(predictions == -1, 1, 0)
        elif method == "lof":
            model.fit(X_scaled)
            scores = model.negative_outlier_factor_
            threshold = np.percentile(scores, 100 * self.contamination)
            predictions = np.where(scores < threshold, 1, 0)
        elif method == "dbscan":
            model.fit(X_scaled)
            predictions = np.where(model.labels_ == -1, 1, 0)
        else:
            predictions = np.zeros(X.shape[0])

        return predictions

    def plotly_outliers(self, X, predictions, sample_names):
        """Plots the outliers using Plotly for interactive visualization.

        Args:
            X (numpy.ndarray): The data points.
            predictions (numpy.ndarray): The outlier predictions.
            sample_names (list): List of sample names.
        """
        X_reduced = PCA(n_components=2).fit_transform(X)

        fig = go.Figure()

        # Plot inliers
        inliers = X_reduced[predictions == 0]
        fig.add_trace(
            go.Scatter(
                x=inliers[:, 0],
                y=inliers[:, 1],
                mode="markers",
                name="Inliers",
                marker=dict(color="blue"),
            )
        )

        # Plot outliers
        outliers = X_reduced[predictions == 1]
        fig.add_trace(
            go.Scatter(
                x=outliers[:, 0],
                y=outliers[:, 1],
                mode="markers",
                name="Outliers",
                marker=dict(color="red"),
            )
        )

        # Add annotations for each sample
        for i, txt in enumerate(sample_names):
            fig.add_annotation(
                x=X_reduced[i, 0], y=X_reduced[i, 1], text=txt, showarrow=False
            )

        # Update layout for better visualization
        fig.update_layout(
            title="Outlier Detection with PCA",
            xaxis_title="Principal Component 1",
            yaxis_title="Principal Component 2",
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        )

        fig.show()

    def reduce_dimensionality(self, X):
        pca = PCA(n_components=2)
        X_reduced = pca.fit_transform(X)
        return X_reduced

    def plot_outliers(self, X, predictions, sample_names):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        # Plot before removing outliers
        X_reduced = self.reduce_dimensionality(X)
        ax[0].scatter(X_reduced[:, 0], X_reduced[:, 1], label="Inliers", color="b")
        ax[0].scatter(
            X_reduced[predictions == 1][:, 0],
            X_reduced[predictions == 1][:, 1],
            label="Outliers",
            color="r",
        )
        for i, txt in enumerate(sample_names):
            ax[0].annotate(txt, (X_reduced[i, 0], X_reduced[i, 1]))
        ax[0].set_title("Before Outlier Removal")
        ax[0].set_xlabel("Principal Component 1")
        ax[0].set_ylabel("Principal Component 2")
        ax[0].legend()

        # Plot after removing outliers
        X_reduced_no_outliers = X_reduced[predictions == 0]
        ax[1].scatter(
            X_reduced_no_outliers[:, 0],
            X_reduced_no_outliers[:, 1],
            label="Inliers",
            color="g",
        )
        ax[1].set_title("After Outlier Removal")
        ax[1].set_xlabel("Principal Component 1")
        ax[1].set_ylabel("Principal Component 2")

        plt.show()


class EnsembleOutlierDetection(OutlierDetection):
    """
    Enhanced class for outlier detection using ensemble methods and visualization.

    This class inherits from OutlierDetection and adds functionalities for ensemble outlier detection
    and PCA-based visualization.

    Attributes:
        contamination (float): The proportion of outliers in the data.
        scaler (StandardScaler): StandardScaler object for data normalization.

    Examples:
        >>> from staver import EnsembleOutlierDetection
        >>> eod = EnsembleOutlierDetection(contamination=0.05)
        >>> eod.detect_outliers(X)
        >>> eod.outlier_decision(X, threshold=1)
        >>> eod.plot_ensemble_outliers(X, sample_names)
    """

    def __init__(self, contamination=0.05, methods=["iforest", "ocsvm", "lof"]):
        """
        Initializes EnhancedOutlierDetection with a contamination rate.

        Args:
            contamination (float): The proportion of the dataset expected to be outliers.
        """
        super().__init__(contamination)
        self.scaler = StandardScaler()
        self.methods = methods

    def _parallel_outlier_detection(self, X_scaled, method):
        """
        Private helper method for parallel execution of outlier detection algorithms.

        Args:
            X_scaled (numpy.ndarray): The standardized data for outlier detection.
            method (str): The method used for outlier detection.

        Returns:
            numpy.ndarray: Array of outlier predictions (1: outlier, 0: inlier).
        """
        return self.detect_outliers(X_scaled, method)

    def ensemble_outlier_score(self, X):
        """
        Calculates the ensemble outlier score by averaging scores from multiple algorithms.

        Args:
            X (numpy.ndarray): The data for outlier detection.

        Returns:
            numpy.ndarray: The ensemble outlier scores.
        """
        methods = [method.lower() for method in self.methods]
        X_scaled = self.scaler.fit_transform(X)

        # Parallel execution of outlier detection methods
        parallel_results = Parallel(n_jobs=len(methods))(
            delayed(self._parallel_outlier_detection)(X_scaled, method)
            for method in methods
        )

        # Calculate ensemble scores
        ensemble_scores = np.mean(parallel_results, axis=0)
        return ensemble_scores, parallel_results

    def outlier_decision(self, X, threshold=1):
        """
        Determines if a sample is an outlier based on the ensemble outlier score.

        Args:
            X (numpy.ndarray): The data for outlier detection.
            threshold (float): The threshold for classifying a sample as an outlier.

        Returns:
            numpy.ndarray: Decision array where 1 indicates an outlier and 0 indicates a non-outlier.
        """
        ensemble_scores, _ = self.ensemble_outlier_score(X)
        decisions = np.where(ensemble_scores >= threshold, 1, 0)
        return decisions

    def plot_ensemble_outliers(self, X, sample_names=None, threshold=1):
        """
        Plots the ensemble outlier detection results using PCA for visualization.

        Args:
            X (numpy.ndarray): The data points.
            sample_names (list, optional): List of sample names. Defaults to None.
        """
        predictions = self.outlier_decision(X, threshold=threshold)
        pca = PCA(n_components=2)
        X_reduced = pca.fit_transform(self.scaler.fit_transform(X))

        plt.figure(figsize=(6, 5))
        # Ensure predictions are 1D
        predictions = np.squeeze(predictions)
        plt.scatter(
            X_reduced[predictions == 0, 0],
            X_reduced[predictions == 0, 1],
            label="Inliers",
            color="blue",
        )
        plt.scatter(
            X_reduced[predictions == 1, 0],
            X_reduced[predictions == 1, 1],
            label="Outliers",
            color="red",
        )
        if sample_names:
            for i, name in enumerate(sample_names):
                if i < len(X_reduced):
                    plt.text(X_reduced[i, 0], X_reduced[i, 1], name)
        plt.xlabel("Principal Component 1")
        plt.ylabel("Principal Component 2")
        plt.title("Ensemble Outlier Detection with PCA")
        plt.legend()
        plt.show()

    def visualize_peptide_abundance(self, X, sample_names=None):
        """Visualizes the abundance density plots for each sample.

        Args:
            X (numpy.ndarray): The peptide abundance data.
            sample_names (list): List of sample names.
        """
        plt.figure(figsize=(12, 6))
        for i, sample in enumerate(sample_names):
            sns.kdeplot(X[:, i], label=sample)
        plt.title("Peptide Abundance Density Plots")
        plt.xlabel("Abundance")
        plt.ylabel("Density")
        plt.legend()
        plt.show()


def generate_data(contamination=0.00, n_samples=200, n_proteins=7000):
    rng = np.random.RandomState(42)
    X_inliers = rng.normal(
        loc=0, scale=1, size=(int((1 - contamination) * n_samples), n_proteins)
    )
    X_outliers = rng.uniform(
        low=-6, high=6, size=(int(contamination * n_samples), n_proteins)
    )
    X = np.concatenate([X_inliers, X_outliers], axis=0)
    return X


if __name__ == "__main__":
    data = pd.read_csv(
        "/Users/ranpeng/Desktop/github/STAVER/staver/results/CRC_051865_library_abudance.csv",
        index_col=0,
    )
    data.replace(np.nan, data.min().min() / 10, inplace=True)
    data = data.T
    # Create an instance of EnhancedOutlierDetection
    outlier_detector = EnsembleOutlierDetection(
        contamination=0.05, methods=["iforest", "ocsvm", "lof"]
    )

    X_example = data.to_numpy()
    sample_names_example = [f"Sample_{i}" for i in range(1, 11)]

    ensemble_scores, eod = outlier_detector.ensemble_outlier_score(X_example)
    # Perform enhanced outlier detection and visualization
    outlier_detector.plot_ensemble_outliers(
        X_example, sample_names_example, threshold=1
    )
    # outlier_detector.visualize_peptide_abundance(X_example, sample_names_example)
    print(f"Ensemble outlier scores: \n{ensemble_scores}")
    print(f"Outlier decision: {eod}")
