#!/bin/bash
# GSP Local Development Server
# Usage: ./start.sh [port]
# Default port: 8000

set -e

PORT=${1:-8000}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "==> GSP Primer Design Tool"
echo "==> Starting local server on http://127.0.0.1:$PORT"
echo ""

# 1. Make binaries executable
echo "[1/4] Checking binaries..."
chmod +x "$SCRIPT_DIR/bin/"* 2>/dev/null || true
chmod +x "$SCRIPT_DIR/gsp_web/ncbi-blast-2.16.0+/bin/"* 2>/dev/null || true

# 2. Create venv if missing or tied to a different machine's Python
if [ ! -d "$SCRIPT_DIR/venv" ] || ! "$SCRIPT_DIR/venv/bin/python" --version &>/dev/null; then
    echo "[2/4] Creating virtual environment..."
    rm -rf "$SCRIPT_DIR/venv"
    python3 -m venv "$SCRIPT_DIR/venv"
else
    echo "[2/4] Virtual environment exists."
fi

# 3. Install/update dependencies
echo "[3/4] Installing Python dependencies..."
"$SCRIPT_DIR/venv/bin/pip" install -q --upgrade pip
"$SCRIPT_DIR/venv/bin/pip" install -q -r "$SCRIPT_DIR/requirements.txt"

# 4. Run migrations and start server
echo "[4/4] Starting Django server..."
cd "$SCRIPT_DIR/gsp_web"
"$SCRIPT_DIR/venv/bin/python" manage.py migrate --run-syncdb -v 0 2>/dev/null || true
"$SCRIPT_DIR/venv/bin/python" manage.py runserver "127.0.0.1:$PORT"
