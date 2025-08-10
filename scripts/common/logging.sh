#!/usr/bin/env bash

_ts() { date +"%H:%M:%S"; }
log_info()  { echo "[$(hostname)] [$(_ts)] [INFO ] $*"; }
log_warn()  { echo "[$(hostname)] [$(_ts)] [WARN ] $*" >&2; }
log_error() { echo "[$(hostname)] [$(_ts)] [ERROR] $*" >&2; }
die() { log_error "$*"; exit 1; }


