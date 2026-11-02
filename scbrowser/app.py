"""Entry point: panel serve scbrowser/app.py"""
import panel as pn

from scbrowser.ui import build_app

pn.extension("tabulator", sizing_mode="stretch_width", notifications=True)

build_app().servable()
