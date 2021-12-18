Run the file `freq_weighting_observe.m` for the fixed-wing test script. This script only performs a Frequency-domain identification, because each of the 10 experiments has noisy 600,000 time-domain record (causing ill-conditioning and memory issues). 

Read `code_observations.txt` for further notes and details.

Right now, the frequency-domain data is uniformly weighted (W=1) between  6-35 Hz. Adding a non-uniform weighting may increase fitting accuracy. So to may incorporating subsamples of the time-domain chirp signal.

https://www.mathworks.com/help/ident/ug/modal-analysis-of-a-flexible-flying-wing-aircraft.html
