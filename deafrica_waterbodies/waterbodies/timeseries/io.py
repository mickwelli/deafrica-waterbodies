import os
import s3urls
import urllib
import boto3
import botocore
import logging
from botocore.exceptions import ClientError

_log = logging.getLogger(__name__)


def check_if_s3_uri(file_path: str):
    """
    Checks if a file path is an S3 URI.

    Parameters
    ----------
    file_path : str
        File path to check

    Returns
    -------
    bool
        True if the file path is an S3 URI or Object URL.
        False if the file path is a local file path or File URI.
    """

    file_scheme = urllib.parse.urlparse(file_path).scheme

    # Assumption here is urls with the scheme https and http
    # are s3 file paths.

    # All others are assumed to be local files.

    valid_s3_schemes = ["s3", "http", "https"]

    if file_scheme in valid_s3_schemes:
        return True
    else:
        return False


def check_local_dir_exists(
        dir_path: str,
        error_if_exists=True):
    """
    Checks if a specified path is an existing directory.
    if `error_if_exists == True`, raises an error if the directory exists.
    if `error_if_exists == False`, raises an error if the directory does not exist.

    Parameters
    ----------
    dir_path : str
        Path to check.
    error_if_exists: bool, optional
        If True, raise an error if the directory exists.
        If False, raise an error if the directory does NOT exist.
        By default True.
    """
    if error_if_exists:
        if os.path.exists(dir_path):
            if os.path.isdir(dir_path):
                raise FileExistsError(f"Directory {dir_path} exists!")
            else:
                raise NotADirectoryError(f"Directory {dir_path} is not a directory!")
        else:
            pass

    else:
        if os.path.exists(dir_path):
            if not os.path.isdir(dir_path):
                raise NotADirectoryError(f"{dir_path} is not a directory!")
            else:
                pass
        else:
            raise FileNotFoundError(f"Directory {dir_path} does not exist!")


def check_local_file_exists(
        file_path: str,
        error_if_exists=True):
    """
    Checks if a specified path is an existing file.
    if `error_if_exists == True`, raises an error if the file exists.
    if `error_if_exists == False`, raises an error if the file does not exist.

    Parameters
    ----------
    file_path : str
        Path to check.
    error_if_exists : bool, optional
        If True, raise an error if the file exists.
        If False, raise an error if the file does NOT exist.
        By default True.
    """
    if error_if_exists:
        if os.path.exists(file_path):
            if not os.path.isfile(file_path):
                raise ValueError(f"{file_path} is not a file!")
            else:
                raise FileExistsError(f"File {file_path} exists!")
        else:
            pass

    else:
        if os.path.exists(file_path):
            if not os.path.isfile(file_path):
                raise ValueError(f"{file_path} is not a file!")
            else:
                pass
        else:
            raise FileNotFoundError(f"{file_path} does not exist!")
    return 0


def check_s3_bucket_exists(
        bucket_name: str,
        s3_client: botocore.client.S3 = None):
    """
    Check if a bucket exists and if the user has permission to access it.

    Parameters
    ----------
    bucket_name : str
        Name of s3 bucket to check.
    s3_client : botocore.client.S3
        A low-level client representing Amazon Simple Storage Service (S3), by default None.

    """
    # Get the service client.
    if s3_client is None:
        s3_client = boto3.client("s3")

    try:
        response = s3_client.head_bucket(Bucket=bucket_name) # noqa E501
    except ClientError as error:
        error_code = int(error.response['Error']['Code'])

        if error_code == 403:
            raise PermissionError(f"{bucket_name} is a private Bucket. Forbidden Access!")
        elif error_code == 404:
            raise FileNotFoundError(f"Bucket {bucket_name} Does Not Exist!")
    except Exception as error:
        _log.exception(error)
        raise error


def check_s3_object_exists(
        s3_object_uri: str,
        error_if_exists=True,
        s3_client: botocore.client.S3 = None):
    """
    Check if an object in an S3 bucket exists.
    if error_if_exists is True, raises an error if the object exists.
    if error_if_exists is False, raises an error if the object does not exist.

    Parameters
    ----------
    s3_object_uri : str
        S3 URI of the object to check.
    error_if_exists : bool, optional
        If True, raise an error if the object exists.
        If False, raise an error if the object does NOT exist.
        By default True.
    s3_client : botocore.client.S3
        A low-level client representing Amazon Simple Storage Service (S3), by default None.
    """
    # Get the service client.
    if s3_client is None:
        s3_client = boto3.client("s3")

    bucket_name = s3urls.parse_url(s3_object_uri)["bucket"]
    object_key = s3urls.parse_url(s3_object_uri)["key"]

    # First check if bucket exists.
    check_s3_bucket_exists(bucket_name, s3_client)

    if error_if_exists:
        try:
            response = s3_client.head_object(Bucket=bucket_name, Key=object_key)
            raise FileExistsError(f"Object {s3_object_uri} already exists!")
        except ClientError as error:
            error_code = int(error.response['Error']['Code'])

            # Object exists but user does not have access.
            if error_code == 403:
                raise FileExistsError(f"Object {s3_object_uri} already exists! Forbidden Access!")
            # Object does not exist.
            elif error_code == 404:
                pass
        except Exception as error:
            _log.exception(error)
            raise error
    else:
        try:
            response = s3_client.head_object(Bucket=bucket_name, Key=object_key)  # noqa E501
        except ClientError as error:
            error_code = int(error.response['Error']['Code'])

            # Object exists but user does not have access.
            if error_code == 403:
                raise PermissionError(f"Object {s3_object_uri} exists, but forbidden access!")
            # File does not exist.
            elif error_code == 404:
                raise FileNotFoundError(f"Object {s3_object_uri} does not exist!")
        except Exception as error:
            _log.exception(error)
            raise error
