import os
from datetime import datetime

ARCHIVE_ENABLED = config.get("archive", {}).get("enabled", False)
ARCHIVE_FOLDER = config.get("archive", {}).get("folder", "archive")

def get_archive_path():
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    return f"{ARCHIVE_FOLDER}/{timestamp}"