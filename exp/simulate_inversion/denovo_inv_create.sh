#!/bin/bash

set -euo pipefail

wd=$(dirname $0)
asm="${1}"

# https://github.com/lh3/minimap2/issues/106
minimap2 -PD -k19 -w19 -m200 -t8 "${asm}" "${asm}" > "${wd}/out.paf"

