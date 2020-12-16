# MESS
Catalog augmentation (event detection) based on matched filter technique (MFT) <br>
<br>
Usage: <br>
(1) Prepare template phase file & Cut template waveform <br>
(2) Run MESS <br>

# Run MESS  
(1) Matched filter: calculate CC trace on every station's record <br>
(2) Expand peak CC value: expand peak value on CC traces <br>
(3) Shift: time shift to origin times for all stations <br>
(4) Stack & Detect on stacked CC trace <br>
(5) dt_p, dt_s are picked by cross-correlation 
