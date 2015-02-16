set "fN=wam_hb_2"

del %fN%.1 %fN%.2 %fN%.3 %fN%.4 %fN%.5 %fN%.6p %fN%.6vx %fN%.6vy %fN%.6vz %fN%.fpt %fN%.out %fN%.p2f errorf.log errorp.log rgkinit.txt rgklog.txt wamitlog.txt %fN%_batch.log

set "t0=%Time%"
set "d0=%Date%"

N:\wamitv7\wamit

set "t1=%Time%"
set "d1=%Date%"

echo WAMIT Run: %fN% >> %fN%_batch.log
echo Started: %d0% %t0% >> %fN%_batch.log
echo Stopped: %d1% %t1% >> %fN%_batch.log
