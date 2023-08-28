"""
Various I/O adaptors
"""
import os
import boto3
import uuid
import shutil
import logging
from botocore.client import ClientError

_log = logging.getLogger(__name__)


def test_access_to_bucket(bucket_name):

    s3_client = boto3.client('s3')

    try:
        test = s3_client.head_bucket(Bucket="deafrica-waterbodies-dev")
        if test:
            _log.info(f"This user can write to the s3 bucket {bucket_name}")
    except ClientError as error:
        _log.error(f"This user can not write to the s3 bucket {bucket_name}")
        _log.error(error)
        raise
    except Exception as error:
        _log.error(error)
        raise


# From https://boto3.amazonaws.com/v1/documentation/api/latest/guide/s3-uploading-files.html
def upload_file(file_name, bucket_name, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket_name: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = os.path.basename(file_name)

    # Upload the file
    s3_client = boto3.client('s3')
    try:
        response = s3_client.upload_file(file_name, bucket_name, object_name)
    except ClientError as e:
        _log.error(e)
        raise


def write_waterbodies_to_file(
        waterbodies_gdf,
        product_version,
        storage_location,
        output_bucket_name,
        output_local_folder,
        output_file_name,
        output_file_type):

    # Validate output file type.
    valid_output_file_type = ['GeoJSON', 'Shapefile']
    try:
        if output_file_type == "GeoJSON":
            output_file_extension = ".geojson"
        elif output_file_type == "Shapefile":
            output_file_extension = ".shp"
        else:
            raise ValueError(f"{output_file_type} is not implemented. Please select from {valid_output_file_type}.")
    except Exception as error:
        _log.error(error)
        raise

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
            raise

    elif storage_location == "s3":

        test_access_to_bucket(bucket_name=output_bucket_name)

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
                    s3_file_uri = f"s3://{output_bucket_name}/{object_prefix}{local_temp_file}"
                    upload_file(file_name=local_temp_file_fp,
                                bucket_name=output_bucket_name,
                                object_name=s3_file_uri)

                # Delete temporary folder.
                shutil.rmtree(local_temp_dir)

            _log.info(f"Waterbodies written to s3 bucket {output_bucket_name} as {s3_uri}")

        except Exception as error:
            _log.error(error)
            raise
