from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("captnseq-toolkit")
except PackageNotFoundError:  # not installed
    __version__ = "0+unknown"
