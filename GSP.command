#!/bin/bash
# GSP Primer Design Tool — Double-click to launch
# Works on macOS: opens Terminal, starts the server, and opens the browser.

PORT=8000
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "==> GSP Primer Design Tool"
echo "==> Starting on http://127.0.0.1:$PORT"
echo ""

# Make binaries executable
chmod +x "$SCRIPT_DIR/bin/"* 2>/dev/null || true
chmod +x "$SCRIPT_DIR/gsp_web/ncbi-blast-2.16.0+/bin/"* 2>/dev/null || true

# Create venv if missing or if it points to a different machine's Python
if [ ! -d "$SCRIPT_DIR/venv" ] || ! "$SCRIPT_DIR/venv/bin/python" --version &>/dev/null; then
    echo "[1/3] Creating virtual environment..."
    rm -rf "$SCRIPT_DIR/venv"
    python3 -m venv "$SCRIPT_DIR/venv"
else
    echo "[1/3] Virtual environment OK."
fi

# Install/update dependencies
echo "[2/3] Checking dependencies..."
"$SCRIPT_DIR/venv/bin/pip" install -q --upgrade pip
"$SCRIPT_DIR/venv/bin/pip" install -q -r "$SCRIPT_DIR/requirements.txt"

# Run migrations silently
cd "$SCRIPT_DIR/gsp_web"
"$SCRIPT_DIR/venv/bin/python" manage.py migrate --run-syncdb -v 0 2>/dev/null || true

# Open browser after a short delay (gives server time to start)
echo "[3/3] Starting server..."
sleep 2 && open "http://127.0.0.1:$PORT" &

# Start server
"$SCRIPT_DIR/venv/bin/python" manage.py runserver "127.0.0.1:$PORT"
