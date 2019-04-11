import numpy as np
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import RobustScaler, StandardScaler


def load_metrics(file_in):
    """
    Load the metrics for each msrun and sample
    :param file_in: File with the corresponding metrics
    :return: Panda Dataframe Metrics
    """
    metrics = pd.read_csv(file_in, '\t', index_col=0, parse_dates=[1])
    metrics.fillna(0, inplace=True)
    # sort the experiments chronologically
    metrics.sort_values(by='StartTimeStamp', inplace=True)
    metrics.set_index('StartTimeStamp', append=True, inplace=True)
    metrics.columns = [m.replace('_', '-') for m in metrics.columns.values]

    print(metrics.head())

    return metrics


def preprocess(data, min_variance, min_corr, scaling_mode):
    # remove low-variance and correlated metrics
    #data, variance = remove_low_variance_features(data, min_variance)
    variance = 0
    data, corr = remove_correlated_features(data, min_corr)
    # scale the values
    data = scale(data, scaling_mode)

    return data, variance, corr


def remove_low_variance_features(data, min_variance):
    variance_threshold = VarianceThreshold(min_variance).fit(data)

    return data[np.where(variance_threshold.variances_ > min_variance)[0]], variance_threshold.variances_


def remove_correlated_features(data, min_corr):
    corr = data.corr()

    # find strongly correlated features
    remove_features = set()
    for row in range(len(corr.index)):
        if corr.columns.values[row] not in remove_features:
            for col in range(row + 1, len(corr.columns)):
                if corr.columns.values[col] not in remove_features and abs(corr.iloc[row, col]) > min_corr:
                    remove_features.add(corr.columns.values[col])

    # remove the correlated features
    for feature in remove_features:
        del data[feature]

    return data, corr


def scale(data, mode='robust'):
    scaler = StandardScaler() if mode == 'standard' else RobustScaler()
    return pd.DataFrame(scaler.fit_transform(data), index=data.index, columns=data.columns)
