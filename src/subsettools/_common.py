"""Helper functions for the subsettools package."""

from datetime import datetime
import pytz


def get_utc_time(date_string, time_zone):
    """Convert the given date and time_zone to UTC time.

    Args:
        date_string (str): date in the form 'yyyy-mm-dd'
        time_zone (str): a pytz-supported time zone

    Returns:
        A timezone-unaware datetime object representing the time in UTC.
    """
    date = datetime.strptime(date_string, "%Y-%m-%d")
    if time_zone != "UTC":
        date = (
            date.replace(tzinfo=pytz.timezone(time_zone))
            .astimezone(pytz.UTC)
            .replace(tzinfo=None)
        )
    return date
