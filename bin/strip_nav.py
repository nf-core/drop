#!/usr/bin/env python3

# MIT License

# Copyright (c) 2019, Michaela Mueller, Vicente Yepez, Christian Mertes, Daniela Andrade, Julien Gagneur

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Remove <div class="side-nav-wrapper"> sections from a MultiQC HTML file.
"""

import argparse
from pathlib import Path
from bs4 import BeautifulSoup

def main():
    ap = argparse.ArgumentParser(description="Remove <div class='side-nav-wrapper'> from HTML")
    ap.add_argument("--html", required=True, help="Input HTML file")
    ap.add_argument("--output", help="Output HTML file (default: *_clean.html)")
    ap.add_argument("--inplace", action="store_true", help="Overwrite input file")
    args = ap.parse_args()

    in_path = Path(args.html)
    html_text = in_path.read_text(encoding="utf-8", errors="ignore")

    soup = BeautifulSoup(html_text, "lxml")
    for div in soup.select("div.side-nav-wrapper"):
        div.decompose()

    cleaned_html = str(soup)

    if args.inplace:
        in_path.write_text(cleaned_html, encoding="utf-8")
        print(f"✔ Removed side-nav-wrapper. Overwrote {in_path}")
    else:
        out_path = Path(args.output) if args.output else in_path.with_name(in_path.stem + "_clean.html")
        out_path.write_text(cleaned_html, encoding="utf-8")
        print(f"✔ Removed side-nav-wrapper. Wrote {out_path}")

if __name__ == "__main__":
    main()
