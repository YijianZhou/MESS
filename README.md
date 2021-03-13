# MESS
Catalog augmentation (event detection) based on matched filter technique (MFT) <br>
<br>
This method includes four steps: <br>
(1) Match: calculate CC trace on every station (matched filter) <br>
(2) Expand: expand peak values on CC traces <br>
(3) Shift: time shift to origin times for all CC traces <br>
(4) Stack: stack CC traces of different stations & detect events on the stacked trace <br>
(5) dt_p, dt_s are picked by cross-correlation 

Usage (see example_mess_workdir): <br>
(1) Prepare template phase file & Cut template waveform <br>
(2) Run MESS <br>
