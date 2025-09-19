#!/usr/bin/env python3
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
