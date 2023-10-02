#!/bin/bash

# Add the wofs_ls_summary_alltime datasets.
s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x204/y048/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'
s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x204/y049/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'

s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x205/y048/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'
s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x205/y049/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'

s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x206/y048/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'
s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x206/y049/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'