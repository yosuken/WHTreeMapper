#!/bin/bash
set -euo pipefail

repo_dir="$(cd "$(dirname "$0")" && pwd)"
env_file="${repo_dir}/environment.yaml"
env_name="$(awk -F': *' '$1 == "name" { print $2; exit }' "${env_file}")"

if [[ -z "${env_name}" ]]; then
  echo "Error: failed to read environment name from ${env_file}" >&2
  exit 1
fi

if command -v micromamba >/dev/null 2>&1; then
  conda_cmd="micromamba"
elif command -v mamba >/dev/null 2>&1; then
  conda_cmd="mamba"
elif command -v conda >/dev/null 2>&1; then
  conda_cmd="conda"
else
  echo "Error: micromamba, mamba, or conda is required." >&2
  exit 1
fi

if "${conda_cmd}" env list | awk '{ print $1 }' | grep -Fxq "${env_name}"; then
  if [[ "${conda_cmd}" == "micromamba" ]]; then
    "${conda_cmd}" install -n "${env_name}" -f "${env_file}" -y
  else
    "${conda_cmd}" env update -n "${env_name}" -f "${env_file}" -y
  fi
else
  "${conda_cmd}" env create -f "${env_file}" -y
fi

env_prefix="$("${conda_cmd}" run -n "${env_name}" sh -c 'printf "%s" "${CONDA_PREFIX}"')"
launcher="${env_prefix}/bin/wh_tree_mapper"

cat > "${launcher}" <<EOF
#!/bin/sh
exec "${repo_dir}/wh_tree_mapper" "\$@"
EOF
chmod 755 "${launcher}"

echo ""
echo "Done."
echo "Activate with: ${conda_cmd} activate ${env_name}"
echo "Then run: wh_tree_mapper -h"
