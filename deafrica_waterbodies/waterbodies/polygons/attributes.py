def get_timeseries_url_from_uid(uid, timeseries_output_bucket, timeseries_product_version):
    version = timeseries_product_version.replace(".", "-")

    subfolder = uid[:4]

    csv_file = f"{uid}_v{version[0]}.csv"

    timeseries_url = f's3://{timeseries_output_bucket}/{version}/timeseries/{subfolder}/{csv_file}'

    return timeseries_url


def add_attributes(waterbodies_gdf, timeseries_output_bucket, timeseries_product_version):

    # Reproject into a projected crs
    waterbodies_gdf_6933 = waterbodies_gdf.to_crs("EPSG:6933")

    # Perimeter
    waterbodies_gdf_6933["perimeter_m"] = waterbodies_gdf_6933.geometry.length
    # Area
    waterbodies_gdf_6933["area_m2"] = waterbodies_gdf_6933.geometry.area
    # Timeseries url
    waterbodies_gdf_6933["timeseries"] = waterbodies_gdf_6933.apply(lambda row: get_timeseries_url_from_uid(row["UID"], timeseries_output_bucket, timeseries_product_version), axis=1)

    # Project to epsg:4326
    waterbodies_gdf_4326 = waterbodies_gdf_6933.to_crs("EPSG:4326")

    return waterbodies_gdf_4326
