import logging
import urllib
from pathlib import Path

import boto3
import fsspec
from mypy_boto3_s3 import S3Client

_log = logging.getLogger(__name__)


def check_if_s3_uri(file_path: str | Path) -> bool:
    """
    Checks if a file path is an S3 URI.

    Parameters
    ----------
    file_path : str | Path
        File path to check

    Returns
    -------
    bool
        True if the file path is an S3 URI.
    """
    # "Support" pathlib Paths.
    file_path = str(file_path)

    file_scheme = urllib.parse.urlparse(file_path).scheme

    valid_s3_schemes = ["s3"]

    if file_scheme in valid_s3_schemes:
        return True
    else:
        return False


def check_dir_exists(dir_path: str | Path) -> bool:
    """
    Checks if a specified path is an existing directory.

    Parameters
    ----------
    dir_path : str | Path
        Path to check.

    Returns
    -------
    bool
        True if the path exists and is a directory.
        False if the path does not exists or if the path exists and it is not a directory.
    """
    # "Support" pathlib Paths.
    dir_path = str(dir_path)

    if check_if_s3_uri(dir_path):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    if fs.exists(dir_path):
        if fs.isdir(dir_path):
            return True
        else:
            return False
    else:
        return False


def check_file_exists(file_path: str | Path) -> bool:
    """
    Checks if a specified path is an existing file.

    Parameters
    ----------
    file_path : str | Path
        Path to check.

    Returns
    -------
    bool
        True if the path exists and is a file.
        False if the path does not exists or if the path exists and it is not a file.
    """
    # "Support" pathlib Paths.
    file_path = str(file_path)

    if check_if_s3_uri(file_path):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    if fs.exists(file_path):
        if fs.isfile(file_path):
            return True
        else:
            return False
    else:
        return False


def check_s3_bucket_exists(bucket_name: str, s3_client: S3Client = None):
    """
    Check if a bucket exists and if the user has permission to access it.

    Parameters
    ----------
    bucket_name : str
        Name of s3 bucket to check.
    s3_client : S3Client
        A low-level client representing Amazon Simple Storage Service (S3), by default None.

    """
    # Get the service client.
    if s3_client is None:
        s3_client = boto3.client("s3")

    try:
        response = s3_client.head_bucket(Bucket=bucket_name)  # noqa E501
    except ClientError as error:
        error_code = int(error.response["Error"]["Code"])

        if error_code == 403:
            raise PermissionError(f"{bucket_name} is a private Bucket. Forbidden Access!")
        elif error_code == 404:
            raise FileNotFoundError(f"Bucket {bucket_name} Does Not Exist!")
    except Exception as error:
        _log.exception(error)
        raise error
