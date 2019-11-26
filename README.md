# MSMS
Catalog augmentation (event detection) based on matched filter technique  
<br>
Usage: <br>
(1) Cut out template waveform <br>
(2) Run MSMS <br>

## Run MSMS  
(1) Matched filter: calculate CC trace on every station's record <br>
(2) Shift: time shift to origin times for all stations <br>
(3) Mask CC traces: mask triggered times with peak CC values <br>
(4) Stack & Detect on stacked CC trace <br>
