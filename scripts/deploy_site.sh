#!/bin/bash
# Deploy the committed site/ from a git ref (default origin/main) to cryomics.org/GENICEII,
# stamping a fresh ?v= on the script tags so browsers fetch the new JS.
set -euo pipefail
REF=${1:-origin/main}
cd "$(dirname "$0")/.."
git fetch -q origin
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
git archive "$REF" site | tar -x -C "$T"
V=$(git rev-parse --short "$REF")
sed -i -E "s#<script src=\"([a-z]+)\.js(\?v=[0-9a-z]+)?\"></script>#<script src=\"\1.js?v=$V\"></script>#" "$T/site/index.html"
rsync -a --chmod=D755,F644 --exclude '_*' "$T/site/" dreamhost:cryomics.org/GENICEII/
echo "deployed $REF ($V)"
