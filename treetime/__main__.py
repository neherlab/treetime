#!/usr/bin/env python
"""
Stub function and module used as a setuptools entry point.
Based on augur's __main__.py and setup.py
"""

from __future__ import print_function, division, absolute_import
import sys
from treetime import version, make_parser
import matplotlib

matplotlib.use("AGG")

# Entry point for setuptools-installed script and bin/augur dev wrapper.
def main():
    parser = make_parser()

    params = parser.parse_args()
    return_code = params.func(params)

    sys.exit(return_code)


# Run when called as `python -m treetime`, here for good measure.
if __name__ == "__main__":
    main()
