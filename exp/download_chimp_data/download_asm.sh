#/bin/bash

set -euo pipefail

xargs -P 4 -I {} bash -c "wget {}" < chimp_asm.txt
