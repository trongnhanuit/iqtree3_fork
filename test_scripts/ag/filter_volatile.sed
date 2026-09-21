# Drop lines that legitimately differ between two runs of different binaries
# (timestamps, timings, host, binary path, version banner). Used by
# identity.sh on .iqtree and .log files before byte comparison.
/Date and time/d
/CPU time/d
/wall-clock/d
/Command:/d
/Host:/d
/^Time:/d
/Seed:/d
/Kernel:/d
/IQ-TREE version/d
/IQ-TREE multicore version/d
/ second/d
/ secs/d
/ sec)/d
/ sec$/d
/ sec /d
/% CPU/d
/Time: /d
/took .* rounds/d
/^Version /d
/built /d
/analytical-gradients/d
