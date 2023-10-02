#!/bin/bash

# Add the wofs_ls_summary_alltime datasets.
s3-to-dc "s3://deafrica-services/wofs_ls_summary_alltime/1-0-0/x164/y098/*/*.json" --stac --no-sign-request --skip-lineage 'wofs_ls_summary_alltime'
