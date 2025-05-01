"""Helper functions for the subsettools package."""

from datetime import datetime
from zoneinfo import ZoneInfo
import importlib
import hf_hydrodata as hf


def get_hf_gridded_data(options):
    "Wrapper around hf_hydrodata.get_gridded_data to handle various exceptions."

    try:
        _add_subsettool_version(options)
        data = hf.get_gridded_data(options)
    except ValueError as err:
        raise ValueError(
            f"No HydroData entry found for the requested filters: {options}"
        ) from err
    except Exception as exc:
        raise RuntimeError(
            f"Failed to get data for the requested filters: {options}"
        ) from exc
    return data


def get_utc_time(date_string, time_zone):
    """Convert the given date and time_zone to UTC time.

    Args:
        date_string (str): date in the form 'yyyy-mm-dd'
        time_zone (str): a zoneinfo-supported time zone

    Returns:
        A timezone-unaware datetime object representing the time in UTC.
    """
    date = datetime.strptime(date_string, "%Y-%m-%d")
    if time_zone != "UTC":
        date = (
            date.replace(tzinfo=ZoneInfo(time_zone))
            .astimezone(ZoneInfo('UTC'))
            .replace(tzinfo=None)
        )
    return date

SUBSETTOOLS_VERSION = None
def _add_subsettool_version(options):
    """Add subsettools version to options"""
    
    global SUBSETTOOLS_VERSION
    if SUBSETTOOLS_VERSION is None:
        SUBSETTOOLS_VERSION = importlib.metadata.version("subsettools")
    options["subsettools"] = SUBSETTOOLS_VERSION