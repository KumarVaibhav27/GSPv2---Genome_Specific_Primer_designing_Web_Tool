# GSP v2.0 — Avena Server Deployment Guide

**App:** GSP (Genome Specific Primers) v2.0 — Django web application  
**Server:** avena.pw.usda.gov (Linux)  
**GitHub repo:** https://github.com/KumarVaibhav27/GSPv2---Genome_Specific_Primer_designing_Web_Tool  
**Contact:** Kumar Vaibhav (kvaibhav@umd.edu)

---

## What This App Does

GSP is a Python/Django web application for designing genome-specific PCR and KASP primers in polyploid wheat species. It uses BLAST+, MUSCLE, and Primer3 internally.

---

## Prerequisites Already on Avena

The following are expected to already be present on the server:

| Tool | Expected path |
|---|---|
| `blastn` | `/usr/bin/blastn` |
| `blastdbcmd` | `/usr/bin/blastdbcmd` |
| MUSCLE v5 binary | `/great/home/pksharma/muscle-linux-x86.v5.3` |
| Primer3 binary | `/great/home/pksharma/primer3/src/primer3_core` |
| BLAST genome databases | `/great/bioshare/genomes/` |
| GFF annotation files | `/great/bioshare/genomes/gff/` |

> **Important:** MUSCLE and Primer3 are in pksharma's home directory. Please verify these paths are accessible from the deployment account (readable + executable). If not, copy them to a shared location (e.g. `/great/bioshare/tools/`) and note the new paths — they will be needed in Step 4.

---

## Step 1 — Check Python version

GSP requires Python 3.11 or newer.

```bash
python3 --version
```

If Python 3.11+ is not available, install it or load the appropriate module:

```bash
# Example on HPC-style servers (use whichever is available):
module load python/3.11
```

---

## Step 2 — Clone the repository

Choose an appropriate installation directory (e.g. `/great/home/<your_username>/` or a shared apps folder):

```bash
cd /great/home/<your_username>
git clone https://github.com/KumarVaibhav27/GSPv2---Genome_Specific_Primer_designing_Web_Tool GSP
cd GSP
```

---

## Step 3 — Create a Python virtual environment and install packages

Only 3 Python packages are required:

```bash
python3 -m venv venv
source venv/bin/activate
pip install django openpyxl anthropic
```

Verify the install:

```bash
python -c "import django; print(django.__version__)"
```

---

## Step 4 — Create the environment file

The app reads configuration from `gsp_web/.env`. Create this file with the following content:

```bash
cat > gsp_web/.env << 'EOF'
ANTHROPIC_API_KEY=<provided_by_vaibhav>
DJANGO_SECRET_KEY=<generate_a_random_50_char_string>
DJANGO_DEBUG=false
DJANGO_ALLOWED_HOSTS=avena.pw.usda.gov,10.200.55.17,localhost
EOF
```

**Generate a secret key** (run this and paste the output as DJANGO_SECRET_KEY):

```bash
python3 -c "import secrets; print(secrets.token_urlsafe(50))"
```

**If MUSCLE or Primer3 paths differ from the defaults**, add these lines too:

```
GSP_MUSCLE_PATH=/path/to/muscle_binary
GSP_PRIMER3_PATH=/path/to/primer3_core
```

**If the BLAST database root is different from `/great/bioshare/genomes/`**, add:

```
GSP_BLAST_DB_ROOT=/actual/path/to/blast_databases
```

---

## Step 5 — Verify tool paths

Before starting the server, confirm all binaries exist and are executable:

```bash
# Should each print the tool's version or usage, not "command not found"
blastn -version
blastdbcmd -version

# Read the MUSCLE and Primer3 paths from settings (or use your custom paths):
/great/home/pksharma/muscle-linux-x86.v5.3 --version
/great/home/pksharma/primer3/src/primer3_core --version

# Confirm genome databases are present
ls /great/bioshare/genomes/ | head -5
```

---

## Step 6 — Initialize the Django database

```bash
source venv/bin/activate
cd gsp_web
python manage.py migrate
```

Expected output: several lines of "OK" messages. No errors should appear.

---

## Step 7 — Start the server

### Option A — Quick test (foreground, port 8000)

```bash
source venv/bin/activate
cd /great/home/<your_username>/GSP/gsp_web
python manage.py runserver 0.0.0.0:8000
```

Open a browser and go to: `http://avena.pw.usda.gov:8000`

The app should load with three tabs: "Paste your own sequences", "BLAST against genome", and "KASP markers".

Press `Ctrl+C` to stop.

---

### Option B — Background process (persistent, survives logout)

```bash
source venv/bin/activate
cd /great/home/<your_username>/GSP/gsp_web
nohup python manage.py runserver 0.0.0.0:8000 > ~/gsp_server.log 2>&1 &
echo "GSP server PID: $!"
```

To check if it is running:

```bash
ps aux | grep manage.py
```

To stop it:

```bash
kill <PID>
```

---

### Option C — systemd service (recommended for production, keeps running after reboot)

Create the service file (requires sudo):

```bash
sudo nano /etc/systemd/system/gsp.service
```

Paste the following (replace `<your_username>` and the GSP path):

```ini
[Unit]
Description=GSP v2.0 Django Web App
After=network.target

[Service]
User=<your_username>
WorkingDirectory=/great/home/<your_username>/GSP/gsp_web
ExecStart=/great/home/<your_username>/GSP/venv/bin/python manage.py runserver 0.0.0.0:8000
Restart=on-failure
RestartSec=5
EnvironmentFile=/great/home/<your_username>/GSP/gsp_web/.env

[Install]
WantedBy=multi-user.target
```

Enable and start:

```bash
sudo systemctl daemon-reload
sudo systemctl enable gsp
sudo systemctl start gsp
sudo systemctl status gsp
```

---

## Step 8 — Firewall / port access

If the server has a firewall, port 8000 must be open for TCP:

```bash
sudo firewall-cmd --permanent --add-port=8000/tcp
sudo firewall-cmd --reload
```

Or via iptables:

```bash
sudo iptables -A INPUT -p tcp --dport 8000 -j ACCEPT
```

---

## Folder structure after deployment

```
GSP/
├── venv/                        ← Python virtual environment
├── gsp_python/                  ← Core pipeline code (no changes needed)
├── gsp_web/
│   ├── .env                     ← Created in Step 4 (contains API key)
│   ├── manage.py
│   ├── db.sqlite3               ← Created after migrate (Step 6)
│   ├── media/jobs/              ← Auto-created when jobs run
│   └── primer_design/
│       └── templates/
└── bin/                         ← MUSCLE/Primer3 binaries (macOS versions)
                                    Not used on Linux — paths set via .env
```

---

## Quick troubleshooting

| Problem | Likely cause | Fix |
|---|---|---|
| `ImportError: No module named django` | venv not activated | Run `source venv/bin/activate` first |
| `DisallowedHost` error in browser | Server IP not in ALLOWED_HOSTS | Add the IP to `DJANGO_ALLOWED_HOSTS` in `.env` |
| BLAST jobs fail / no results | Wrong `GSP_BLAST_DB_ROOT` | Check the path in `.env` matches actual DB location |
| MUSCLE not found | Path in `.env` is wrong or not executable | Run `chmod +x <muscle_path>` and verify path |
| Primer3 not found | Same as above | Run `chmod +x <primer3_path>` and verify path |
| Port 8000 not reachable | Firewall blocking | See Step 8 |
| AI chatbot not responding | Missing or wrong `ANTHROPIC_API_KEY` | Verify key in `.env` |

---

## Summary — minimum required actions

1. `git clone` the repo
2. `pip install django openpyxl anthropic` inside a venv
3. Create `gsp_web/.env` with the API key and secret key
4. `python manage.py migrate`
5. `python manage.py runserver 0.0.0.0:8000`

Everything else (BLAST databases, binaries) is already on Avena.
