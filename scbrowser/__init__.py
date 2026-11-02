"""Single-cell gene browser (Panel + Bokeh)."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scbrowser")
except PackageNotFoundError:   # running from a source tree that was not pip-installed
    __version__ = "0.1.0"
