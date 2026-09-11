from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("AAFTF")
except PackageNotFoundError:
    # Defensive fallback in case someone imports directly from source without installing
    __version__ = "0.0.0+unknown"
