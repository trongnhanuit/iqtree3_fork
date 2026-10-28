#!/usr/bin/env python3
"""Bake an --ag-trace file into one self-contained HTML page that opens offline.

usage: bake_trace.py TRACE[.gz] [-o OUT.html] [--lib plotly.min.js]

Plotly is not kept in git: it is downloaded once (pinned version, SHA-256
checked) and cached in ~/.cache/iqtree-agviz/ (or $XDG_CACHE_HOME). --lib uses
a local copy instead (still checked).
"""
import argparse
import base64
import gzip
import hashlib
import os
import subprocess
import sys
import urllib.request

URL = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"
SHA256 = "6d21266ce1bd7d9e5ab4e115989c70c20de0382fd973a8f26ab58619eba4d603"
HERE = os.path.dirname(os.path.abspath(__file__))


def library(path):
    if not path:
        cache = os.path.join(os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "iqtree-agviz")
        path = os.path.join(cache, "plotly-2.35.2.min.js")
        if not os.path.exists(path):
            os.makedirs(cache, exist_ok=True)
            print("downloading " + URL, file=sys.stderr)
            try:
                data = urllib.request.urlopen(URL).read()
                with open(path + ".tmp", "wb") as f:
                    f.write(data)
            except Exception:  # e.g. CERTIFICATE_VERIFY_FAILED on python.org builds: try curl
                if subprocess.call(["curl", "-sSfL", "-o", path + ".tmp", URL]) != 0:
                    sys.exit("download failed; fetch %s yourself and pass it with --lib" % URL)
            os.replace(path + ".tmp", path)
    data = open(path, "rb").read()
    if hashlib.sha256(data).hexdigest() != SHA256:
        sys.exit(path + ": SHA-256 does not match the pinned Plotly 2.35.2; refusing to inline it")
    return data.decode("utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("trace")
    ap.add_argument("-o", dest="out")
    ap.add_argument("--lib")
    args = ap.parse_args()
    trace = open(args.trace, "rb").read()
    if not trace.startswith(b"\x1f\x8b"):
        trace = gzip.compress(trace, mtime=0)
    html = open(os.path.join(HERE, "visualizer.html"), encoding="utf-8").read()
    start, end = "<!-- AGVIZ:LIB -->", "<!-- /AGVIZ:LIB -->"
    i, j = html.index(start), html.index(end) + len(end)
    html = html[:i] + "<script>" + library(args.lib) + "</script>" + html[j:]
    data = '<script id="trace-data" type="application/octet-stream">' + base64.b64encode(trace).decode() + "</script>"
    html = html.replace("<!-- AGVIZ:DATA -->", data)
    out = args.out or args.trace.replace(".gz", "").replace(".jsonl", "") + ".html"
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print("wrote " + out)


if __name__ == "__main__":
    main()
