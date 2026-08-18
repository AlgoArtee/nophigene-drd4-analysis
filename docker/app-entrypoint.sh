#!/bin/sh
set -eu

if [ -f "${NOPHIGENE_DATABASE_PATH}" ]; then
  python src/app.py db backup --label daily --password-file "${NOPHIGENE_DATABASE_KEY_FILE}"
  if [ "$(date +%u)" = "7" ]; then
    python src/app.py db backup --label weekly --password-file "${NOPHIGENE_DATABASE_KEY_FILE}"
  fi
fi

alembic upgrade head
exec python src/app.py "$@"
