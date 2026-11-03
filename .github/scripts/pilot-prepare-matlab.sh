#!/usr/bin/env bash
set -euo pipefail

# Bound dependency downloads before setup-matlab invokes apt. Some hosted
# runners advertise an unreachable Azure HTTP mirror and otherwise retry
# package indexes for most of the job. Keep Ubuntu's signed repositories.
if [[ "${RUNNER_OS:-Linux}" != Linux ]]; then exit 0; fi
if [[ -f /etc/apt/apt-mirrors.txt ]]; then
  printf '%s\n' 'https://archive.ubuntu.com/ubuntu' | sudo tee /etc/apt/apt-mirrors.txt >/dev/null
fi
printf '%s\n' \
  'Acquire::Retries "2";' \
  'Acquire::http::Timeout "20";' \
  'Acquire::https::Timeout "20";' \
  | sudo tee /etc/apt/apt.conf.d/99-pipecompare-timeouts >/dev/null
echo 'MATLAB dependency downloads use bounded timeouts.'
