set sysroot /
set mi-async off
target extended-remote 127.0.0.1:2650
source /home/bolker/Documents/R/pkgs/lme4/misc/asan_rr_debug/rr_watch.py
continue
