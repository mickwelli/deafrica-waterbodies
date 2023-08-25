import os
import math
import logging
import click
import shapely
import geopandas as gpd

from .main import main, logging_setup

from deafrica_waterbodies.waterbodies.polygons.make_polygons import get_waterbodies


@main.command("waterbodies-from-boundingbox", 
              short_help="Waterbodies for area defined by bounding box.",
              no_args_is_help=True)
@click.option("--bbox",
              help="Coordinates of the area of interest's bounding box in the format 'minx,miny,maxx,maxy'.")
@click.option("--bbox-crs",
              default="EPSG:4326",
              help="CRS of the bounding box coordinates.",
              show_default=True)
@click.option("--primary-threshold",
              default=0.1,
              type=click.FLOAT,
              help="Threshold to the location of the waterbody polygons.",
              show_default=True)
@click.option("--secondary-threshold",
              default=0.05,
              type=click.FLOAT,
              help="Threshold to define the extent/shape of the waterbody polygons.",
              show_default=True)
@click.option("--min-polygon-size",
              default=4500,
              show_default=True,
              type=click.FLOAT,
              help="Minimum area in m2 of the waterbody polygons to be included.")
@click.option("--max-polygon-size",
              type=click.FLOAT,
              default=math.inf,
              show_default=True,
              help="Maximum area in m2 of the waterbody polygons to be included.")
@click.option("--remove-ocean-polygons",
              is_flag=True,
              default=False,
              show_default=True,
              help="Filter out ocean polygons from the waterbody polygons.")
@click.option("--land-sea-mask-fp",
              type=click.Path(),
              default=None,
              help="File path to vector dataset to use to filter out ocean polygons.")
@click.option("--remove-major-rivers",
              is_flag=True,
              default=False,
              help="Filter out major rivers from the waterbody polygons.")
@click.option("--major-rivers-mask-fp",
              type=click.Path(),
              default=None,
              help="File path to vector dataset to use to filter out major rivers.")
@click.option("--remove-cbd",
              is_flag=True,
              default=False,
              help="Filter out CBDs from the waterbody polygons.")
@click.option("--urban-mask-fp",
              type=click.Path(),
              default=None,
              help="File path to vector dataset to use to filter out urban/CBD areas.")
@click.option("--handle-large-polygons",
              default="nothing",
              type=click.Choice(['erode-dilate-v1', 'erode-dilate-v2', 'nothing']),
              show_default=True,
              help="Method to use to split large polygons above the Polsby–Popper test value."
              )
@click.option("--pp-test-threshold",
              type=click.FLOAT,
              default=0.005,
              show_default=True,
              help="Polsby–Popper test value to use to split large polygons.")
@click.option("-v", "--verbose", count=True)
@click.option("--product-version",
              type=str,
              default="0.0.1",
              show_default=True,
              help="Product version for the DE Africa Waterbodies product.")
@click.option("--s3/--local",
              default=False,
              help="Save the output to an s3 bucket or a local folder.")
@click.option("--ouptut-bucket-name",
              type=str,
              default="deafrica-waterbodies-dev",
              show_default=True,
              required=False,
              is_eager=True,
              help="The s3 bucket to write the output to.",)
@click.option("--ouptut-local-folder",
              type=click.Path(),
              default=os.getcwd(),
              show_default="Current working directory.",
              required=False,
              is_eager=True,
              help="Directory to write the waterbody polygons to.",)
@click.option("--output-file-name",
              default="waterbodies",
              show_default=True,
              help="File name for the output waterbodies.")
@click.option("--output-file-type",
              default='GeoJSON',
              show_default=True,
              help="File type for the outputs waterbodies.",
              type=click.Choice(['GeoJSON', 'Shapefile'], case_sensitive=False))
def waterbodies_from_bbox(bbox,
                          bbox_crs,
                          secondary_threshold,
                          primary_threshold,
                          min_polygon_size,
                          max_polygon_size,
                          remove_ocean_polygons,
                          land_sea_mask_fp,
                          remove_major_rivers,
                          major_rivers_mask_fp,
                          remove_cbd,
                          urban_mask_fp,
                          handle_large_polygons,
                          pp_test_threshold,
                          verbose,
                          product_version,
                          s3,
                          output_bucket_name,
                          output_local_folder,
                          output_file_name,
                          output_file_type,
                          ):
    """
    Generate waterbodies for WOfS All Time Summary regions covering the area defined in the provided bounding box.
    """
    logging_setup(verbose)
    _log = logging.getLogger(__name__)

    output_fp, log_msg = setup_output_fp(
        product_version,
        s3,
        output_bucket_name,
        output_local_folder,
        output_file_name,
        output_file_type)
    _log.info(log_msg)

    if remove_ocean_polygons:
        filter_out_ocean_polygons = True
    else:
        filter_out_ocean_polygons = False

    if remove_major_rivers:
        filter_out_major_rivers_polygons = True
    else:
        filter_out_major_rivers_polygons = False

    if remove_cbd:
        filter_out_urban_polygons = True
    else:
        filter_out_urban_polygons = False

    dask_chunks = {"x": 3000, "y": 3000, "time": 1}
    resolution = (-30, 30)
    output_crs = "EPSG:6933"
    min_valid_observations = 128

    # Convert bounding box to GeoDataFrame.
    bbox_ = [float(i.strip()) for i in bbox.split(",")]
    aoi_gdf = gpd.GeoDataFrame(geometry=[shapely.geometry.box(*bbox_)], crs=bbox_crs)
    aoi_gdf = aoi_gdf.to_crs(output_crs)
    continental_run = False

    waterbodies_gdf = get_waterbodies(
        aoi_gdf=aoi_gdf,
        continental_run=continental_run,
        dask_chunks=dask_chunks,
        resolution=resolution,
        output_crs=output_crs,
        min_valid_observations=min_valid_observations,
        primary_threshold=primary_threshold,
        secondary_threshold=secondary_threshold,
        min_polygon_size=min_polygon_size,
        max_polygon_size=max_polygon_size,
        filter_out_ocean_polygons=filter_out_ocean_polygons,
        land_sea_mask_fp=land_sea_mask_fp,
        filter_out_major_rivers_polygons=filter_out_major_rivers_polygons,
        major_rivers_mask_fp=major_rivers_mask_fp,
        filter_out_urban_polygons=filter_out_urban_polygons,
        urban_mask_fp=urban_mask_fp,
        handle_large_polygons=handle_large_polygons,
        pp_test_threshold=pp_test_threshold,
        )

    # Reproject to EPSG:4326
    waterbodies_gdf_4326 = waterbodies_gdf.to_crs("EPSG:4326")

    _log.info(f"Writing waterbodies to {output_fp} ...")
    try:
        # If writing to s3 bucket, 
        # if user has write access to bucket this should work.
        waterbodies_gdf_4326.to_file(output_fp)
        _log.info("Done.")
    except Exception as error:
        _log.error(error)
        raise
