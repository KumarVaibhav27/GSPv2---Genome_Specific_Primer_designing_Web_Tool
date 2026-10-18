#!/bin/bash
# GSP Production Server
# Usage: ./start_server.sh [port] [workers]
# Default port: 8000, workers: 4
#
# Requires gunicorn (installed automatically).
# Set environment variables in .env or export them before running.
# Example:
#   export DJANGO_ALLOWED_HOSTS="yourserver.example.com"
#   export DJANGO_SECRET_KEY="your-secret-key-here"
#   ./start_server.sh 8000 4

set -e

PORT=${1:-8000}
WORKERS=${2:-4}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "==> GSP Primer Design Tool — Production Server"
echo "==> Binding to 0.0.0.0:$PORT with $WORKERS workers"
echo ""

# 1. Make binaries executable
echo "[1/5] Checking binaries..."
chmod +x "$SCRIPT_DIR/bin/"* 2>/dev/null || true
chmod +x "$SCRIPT_DIR/gsp_web/ncbi-blast-2.16.0+/bin/"* 2>/dev/null || true

# 2. Create venv if needed
if [ ! -d "$SCRIPT_DIR/venv" ]; then
    echo "[2/5] Creating virtual environment..."
    python3 -m venv "$SCRIPT_DIR/venv"
else
    echo "[2/5] Virtual environment exists."
fi

# 3. Install dependencies + gunicorn
echo "[3/5] Installing Python dependencies..."
"$SCRIPT_DIR/venv/bin/pip" install -q --upgrade pip
"$SCRIPT_DIR/venv/bin/pip" install -q -r "$SCRIPT_DIR/requirements.txt"
"$SCRIPT_DIR/venv/bin/pip" install -q gunicorn

# 4. Collect static files and run migrations
echo "[4/5] Preparing Django..."
cd "$SCRIPT_DIR/gsp_web"
DJANGO_DEBUG=false "$SCRIPT_DIR/venv/bin/python" manage.py migrate --run-syncdb -v 0 2>/dev/null || true
DJANGO_DEBUG=false "$SCRIPT_DIR/venv/bin/python" manage.py collectstatic --noinput -v 0 2>/dev/null || true

# 5. Start gunicorn
echo "[5/5] Starting gunicorn..."
DJANGO_DEBUG=false "$SCRIPT_DIR/venv/bin/gunicorn" \
    --workers "$WORKERS" \
    --bind "0.0.0.0:$PORT" \
    --timeout 300 \
    --access-logfile - \
    gsp_web.wsgi:application
