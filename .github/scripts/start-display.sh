#!/usr/bin/env bash
set -euo pipefail

# Hosted runners normally have Xvfb already. Do not refresh package indexes
# unless installation is needed: an unavailable mirror must not consume the job.
if ! command -v Xvfb >/dev/null 2>&1; then
  apt_options=(-o Acquire::Retries=2 -o Acquire::http::Timeout=20 -o Acquire::https::Timeout=20)
  sudo timeout 240s apt-get "${apt_options[@]}" update -q
  sudo timeout 240s apt-get "${apt_options[@]}" install -y -q --no-install-recommends xvfb
fi

pc_display="${DISPLAY:-:99}"
pc_display_log="${RUNNER_TEMP:-/tmp}/pipecompare-xvfb.log"
Xvfb "$pc_display" -screen 0 1920x1080x24 >"$pc_display_log" 2>&1 &
pc_xvfb_pid=$!

for ((attempt = 0; attempt < 15; attempt++)); do
  if ! kill -0 "$pc_xvfb_pid" 2>/dev/null; then
    cat "$pc_display_log"
    echo 'Xvfb exited before the display was ready.' >&2
    exit 1
  fi
  pc_ready=false
  if command -v xdpyinfo >/dev/null 2>&1; then
    if timeout 1s xdpyinfo -display "$pc_display" >/dev/null 2>&1; then pc_ready=true; fi
  elif [[ -S "/tmp/.X11-unix/X${pc_display#:}" ]]; then
    pc_ready=true
  fi
  if "$pc_ready"; then
    echo "DISPLAY=$pc_display" >>"$GITHUB_ENV"
    echo "Virtual display $pc_display is ready."
    exit 0
  fi
  sleep 1
done

cat "$pc_display_log"
echo 'Xvfb was not ready after 15 checks.' >&2
kill "$pc_xvfb_pid" 2>/dev/null || true
exit 1
