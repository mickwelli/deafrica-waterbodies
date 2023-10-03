import logging
import os
import shutil
import urllib
import uuid
from pathlib import Path

import boto3
import fsspec
import geopandas as gpd
from botocore.client import ClientError
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


# From https://boto3.amazonaws.com/v1/documentation/api/latest/guide/s3-uploading-files.html
def upload_file_to_s3(
    file_name: str, bucket_name: str, object_name: str = None, s3_client: S3Client = None
) -> bool:
    """
    Upload a file to an S3 bucket.


    Parameters
    ----------
    file_name : str
        File path of the file to upload.
    bucket_name : str
        Name of the s3 bucket to upload to.
    object_name : str, optional
        S3 object name. If not specified then `file_name` is used, by default None
    s3_client : S3Client
        A low-level client representing Amazon Simple Storage Service (S3), by default None.

    Returns
    -------
    bool
        True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = os.path.basename(file_name)

    # Get the service client.
    if s3_client is None:
        s3_client = boto3.client("s3")

    # Upload the file
    try:
        response = s3_client.upload_file(file_name, bucket_name, object_name)  # noqa F841
    except ClientError as error:
        _log.exception(error)
        return False
    else:
        return True


def write_waterbodies_to_file(
    waterbodies_gdf: gpd.GeoDataFrame,
    product_version: str,
    storage_location: str,
    output_bucket_name: str,
    output_local_folder: str,
    output_file_name: str,
    output_file_type: str,
):
    """
    Function to GeoDataFrame of waterbody polygons to an ESRI Shapefile or GeoJSON file.

    Parameters
    ----------
    waterbodies_gdf : gpd.GeoDataFrame
        The waterbody polygons.
    product_version : str
        The DE Africa Waterbodies product version number.
    storage_location : str
        Type of storage location. Either "local" or "s3".
    output_bucket_name : str
        The name of the S3 Bucket to write the polygons to.
    output_local_folder : str
        The file path of the local folder to write the polygons to.
    output_file_name : str
        The name of the file to write the waterbody polygons to.
    output_file_type : str
        File type to write the polygons to. Either "GeoJSON" or "ESRI Shapefile"

    """

    # Validate output file type.
    valid_output_file_type = ["GeoJSON", "Shapefile"]
    try:
        if output_file_type == "GeoJSON":
            output_file_extension = ".geojson"
        elif output_file_type == "Shapefile":
            output_file_extension = ".shp"
        else:
            _log.error(f"{output_file_type} is not implemented.")
            raise ValueError(
                f"Invalid output file type. Select a valid output from {valid_output_file_type}."
            )
    except Exception as error:
        _log.exception(error)
        raise error

    object_prefix = f'{product_version.replace(".", "-")}/shapefile/'
    object_name = f"{output_file_name}{output_file_extension}"

    if storage_location == "local":
        try:
            local_dir_fp = os.path.join(output_local_folder, object_prefix)

            if not os.path.exists(local_dir_fp):
                os.makedirs(local_dir_fp)
                _log.info(f"Output folder {local_dir_fp} created.")

            local_file_fp = os.path.join(local_dir_fp, object_name)
            absolute_local_file_fp = os.path.abspath(local_file_fp)

            waterbodies_gdf.to_file(absolute_local_file_fp)

            _log.info(f"Waterbodies written to local disk as {absolute_local_file_fp}")

        except Exception as error:
            _log.error(error)
            raise error

    elif storage_location == "s3":
        # Get the service client.
        s3_client = boto3.client("s3")

        check_s3_bucket_exists(output_bucket_name, s3_client)

        s3_uri = f"s3://{output_bucket_name}/{object_prefix}{object_name}"

        try:
            if output_file_extension == ".geojson":
                waterbodies_gdf.to_file(s3_uri)

            elif output_file_extension == ".shp":
                # Make a temporary folder
                myuuid = str(uuid.uuid1())[:6]
                local_temp_dir = f"temp_{myuuid}"
                os.mkdir(local_temp_dir)

                # Write the waterbodies to the temporary folder.
                local_temp_shapefile_fp = os.path.join(local_temp_dir, object_name)
                waterbodies_gdf.to_file(local_temp_shapefile_fp)

                # Upload each object in the temporary folder into s3.
                local_temp_files = os.listdir(local_temp_dir)
                for local_temp_file in local_temp_files:
                    local_temp_file_fp = os.path.join(local_temp_dir, local_temp_file)
                    upload_file_to_s3(
                        file_name=local_temp_file_fp,
                        bucket_name=output_bucket_name,
                        object_name=f"{object_prefix}{local_temp_file}",
                        s3_client=s3_client,
                    )

                # Delete temporary folder.
                shutil.rmtree(local_temp_dir)

            _log.info(f"Waterbodies written to s3 bucket {output_bucket_name} as {s3_uri}")

        except Exception as error:
            _log.error(error)
            raise error
