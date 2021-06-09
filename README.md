# MESS
MESS is a matched filter earthquake detector <br>
<br>
This method includes four steps: <br>
(1) **M**atch: calculate CC trace on every station (matched filter) <br>
(2) **E**xpand: expand peak values on CC traces <br>
(3) **S**hift: time shift to origin times for all CC traces <br>
(4) **S**tack: stack CC traces of different stations & detect events on the stacked trace <br>
(5) *dt_p* and *dt_s* are picked by cross-correlation 

Usage (see __example_mess_workdir__): <br>
(1) Prepare template phase file & Cut template waveform <br>
(2) Run MESS <br>

**Installation** <br>
MESS is a set of codes. All you need is to setup proper Python environment. This can be accomplished easily by installing [Anaconda](https://www.anaconda.com/products/individual#Downloads), [Obspy](https://github.com/obspy/obspy/wiki/Installation-via-Anaconda), and [Pytorch](https://pytorch.org/) sequentially. Or you can use the *env/mess.yml* file in the with *conda*. 
