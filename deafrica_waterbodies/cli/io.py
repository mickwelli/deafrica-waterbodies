"""
Various I/O adaptors
"""
import os
import boto3
import logging

from botocore.client import ClientError

_log = logging.getLogger(__name__)


def test_access_to_bucket(bucket_name):

    s3 = boto3.client('s3')

    try:
        test = s3.head_bucket(Bucket="deafrica-waterbodies-dev")
        if test:
            _log.info(f"This user can write to the s3 bucket {bucket_name}")
    except ClientError as error:
        _log.error(f"This user can not write to the s3 bucket {bucket_name}")
        _log.error(error)
        raise
    except Exception as error:
        _log.error(error)
        raise


def setup_output_fp(
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
            # Output to be saved to local file.
            local_dir_fp = os.path.join(output_local_folder, object_prefix)

            if not os.path.exists(local_dir_fp):
                os.makedirs(local_dir_fp)
                _log.info(f"Output folder {local_dir_fp} created.")

            output_fp = os.path.join(local_dir_fp, object_name)
            _log.info(f"Waterbodies will be written to local disk as {output_fp}")
        except Exception as error:
            _log.error(error)
            raise

    elif storage_location == "s3":

        output_fp = f"s3://{output_bucket_name}/{object_prefix}{object_name}"

        test_access_to_bucket(bucket_name=output_bucket_name)

        _log.info(f"Waterbodies will be written to the s3 bucket {output_bucket_name} as {output_fp}")

    return output_fp
