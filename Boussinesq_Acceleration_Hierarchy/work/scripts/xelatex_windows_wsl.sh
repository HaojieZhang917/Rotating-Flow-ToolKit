#!/usr/bin/env bash
set -euo pipefail

# Run the existing Windows TeX Live from a WSL workspace.  Mapping the UNC
# directory with `pushd` is required because xdvipdfmx cannot use a UNC path as
# its process working directory.
if [[ $# -lt 1 ]]; then
    echo "usage: $0 TEX_FILE" >&2
    exit 2
fi

doc_path="${!#}"
doc_dir="$(cd "$(dirname "$doc_path")" && pwd)"
doc_file="$(basename "$doc_path")"
win_dir="$(wslpath -w "$doc_dir")"
win_xelatex='E:\latex\texlive\2024\bin\windows\xelatex.exe'

/mnt/c/Windows/System32/cmd.exe /d /c \
    "pushd ${win_dir} && ${win_xelatex} -synctex=1 -interaction=nonstopmode -file-line-error -halt-on-error ${doc_file}"
