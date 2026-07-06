#!/usr/bin/env bash
# Run once from the repository root on cmslpc to create the PocketCoffea LPC
# apptainer shell and make this analysis importable inside it.

set -euo pipefail

curl -OL https://raw.githubusercontent.com/PocketCoffea/lpcjobqueue/main/bootstrap.sh
bash bootstrap.sh
rm -f bootstrap.sh

MARKER="# disapptrks-nano LPC setup"

if grep -q "$MARKER" .bashrc 2>/dev/null; then
    echo "DisappTrks_Nano LPC setup already present in .bashrc -- skipping."
else
    cat >> .bashrc << 'EOF'

# disapptrks-nano LPC setup
# The stock pocket-coffea executable in the container can resolve to a Python
# outside the venv. Keep the command on the venv Python that has lpcjobqueue.
if [[ -d .env ]]; then
    printf '#!/usr/bin/env bash\nexec "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/python" -m pocket_coffea "$@"\n' > .env/bin/pocket-coffea
    chmod +x .env/bin/pocket-coffea
fi

# Make the local analysis package importable by workers.
python -m pip install -e '.[analysis]' -q
EOF
    echo "Added DisappTrks_Nano setup to .bashrc."
fi

echo "Done. Run ./shell, then use the LPC Dask command from README.md."
