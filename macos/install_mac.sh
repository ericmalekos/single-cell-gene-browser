#!/usr/bin/env bash
# One-time setup of the gene browser on a Mac, ending with a "Gene Browser" icon on the Desktop.
#
#   bash macos/install_mac.sh [--data /path/to/app_data] [--port 5006]
#
# Run it from a copy of the repository that contains the prebuilt app_data/ folder (or point --data at one).
# Safe to re-run after updating the code or the data: it updates the environment and rebuilds the app.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="gene_browser"
APP_NAME="Gene Browser"
PORT=5006
DATA=""

usage() { sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case "$1" in
    --data) DATA="${2:?--data needs a path}"; shift 2 ;;
    --port) PORT="${2:?--port needs a number}"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown option: $1" >&2; usage 1 ;;
  esac
done
[ "$(uname -s)" = "Darwin" ] || { echo "This installer is for macOS." >&2; exit 1; }
say() { printf '\n==> %s\n' "$*"; }

# 1. conda: an existing installation, or Miniforge installed into ~/miniforge3
CONDA="$(command -v conda || true)"
for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
         "/opt/homebrew/Caskroom/miniforge/base/bin/conda"; do
  [ -n "$CONDA" ] && break
  [ -x "$c" ] && CONDA="$c"
done
if [ -z "$CONDA" ]; then
  say "Installing Miniforge (conda) into ~/miniforge3"
  installer="${TMPDIR:-/tmp}/Miniforge3-installer.sh"
  curl -fL --progress-bar -o "$installer" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"
  bash "$installer" -b -p "$HOME/miniforge3"
  rm -f "$installer"
  CONDA="$HOME/miniforge3/bin/conda"
fi

# 2. the Python environment with the app (no R needed: the Mac uses prebuilt data)
cd "$REPO"
if "$CONDA" env list | awk '{print $1}' | grep -x "$ENV_NAME" > /dev/null; then
  say "Updating the $ENV_NAME environment"
  "$CONDA" env update -n "$ENV_NAME" -f environment.yml
else
  say "Creating the $ENV_NAME environment (takes a few minutes)"
  "$CONDA" env create -n "$ENV_NAME" -f environment.yml
fi
SCBROWSER="$("$CONDA" info --base)/envs/$ENV_NAME/bin/scbrowser"
if [ ! -x "$SCBROWSER" ]; then   # environments kept outside the base prefix
  SCBROWSER="$("$CONDA" env list | awk -v n="$ENV_NAME" '$1 == n {print $NF}')/bin/scbrowser"
fi
[ -x "$SCBROWSER" ] || { echo "scbrowser was not installed into the $ENV_NAME environment" >&2; exit 1; }

# 3. where the prebuilt data is
if [ -n "$DATA" ]; then
  DATA="$(cd "$DATA" && pwd)"
  touch config.local.yaml
  grep -v '^app_data:' config.local.yaml > config.local.yaml.tmp || true
  mv config.local.yaml.tmp config.local.yaml
  printf 'app_data: "%s"\n' "$DATA" >> config.local.yaml
  say "Using the data in $DATA (saved in config.local.yaml)"
fi
status="$("$SCBROWSER" check)"
printf '%s\n' "$status" | sed -n '/^Stores/,/^$/p'
if printf '%s\n' "$status" | grep "not built" > /dev/null; then
  echo "Warning: some datasets are missing; their tabs will say \"not built\". Copy the app_data/ folder built"
  echo "on the Linux machine into $REPO (or rerun with --data /path/to/app_data)."
fi

# 4. the app: a compiled AppleScript launcher with the gene browser icon
say "Creating ~/Applications/$APP_NAME.app"
SUPPORT="$HOME/Library/Application Support/GeneBrowser"
APP="$HOME/Applications/$APP_NAME.app"
mkdir -p "$SUPPORT" "$HOME/Applications" "$HOME/Library/Logs"
sed -e "s|__SCBROWSER__|$SCBROWSER|g" -e "s|__PORT__|$PORT|g" \
    -e "s|__LOG__|$HOME/Library/Logs/GeneBrowser.log|g" -e "s|__PIDFILE__|$SUPPORT/server.pid|g" \
    "$REPO/macos/GeneBrowser.applescript" > "$SUPPORT/launcher.applescript"
rm -rf "$APP"
osacompile -o "$APP" "$SUPPORT/launcher.applescript"

ICONSET="$SUPPORT/GeneBrowser.iconset"
rm -rf "$ICONSET" && mkdir -p "$ICONSET"
for s in 16 32 128 256 512; do
  sips -z "$s" "$s" "$REPO/macos/icon.png" --out "$ICONSET/icon_${s}x${s}.png" > /dev/null
  sips -z $((s * 2)) $((s * 2)) "$REPO/macos/icon.png" --out "$ICONSET/icon_${s}x${s}@2x.png" > /dev/null
done
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/applet.icns"
rm -rf "$ICONSET"
# recent macOS versions also give applets an asset-catalog icon, which would win over applet.icns
rm -f "$APP/Contents/Resources/Assets.car"
/usr/libexec/PlistBuddy -c "Delete :CFBundleIconName" "$APP/Contents/Info.plist" 2> /dev/null || true
codesign --force --deep --sign - "$APP" 2> /dev/null || true   # re-seal after editing the bundle
touch "$APP"

LINK="$HOME/Desktop/$APP_NAME.app"
[ -L "$LINK" ] && rm "$LINK"
[ -e "$LINK" ] || ln -s "$APP" "$LINK"

say "Done. Double-click \"$APP_NAME\" on the Desktop (or in ~/Applications)."
echo "    It opens http://localhost:$PORT/app in your browser. Double-click again to open it again or stop it."
echo "    Log: ~/Library/Logs/GeneBrowser.log"
