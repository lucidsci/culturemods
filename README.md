Oxygen cell culture models.  Currently implements a 1D kinetics model of change in O2 over time in well media from saturation after a step change in oxygen consumption by a cell monolayer at the bottom.

Example commands.


python3 cli.py o2-at-position-by-ocr --position 1.25 --vol 100 50 75 100 125 150 175 200

This will plot the O2 at 1.25mm above thee well bottom for a media volume of 100uL for OCRs of 50 thru 200 fmols/mm2/s

![O2 over time](https://github.com/lucidsci/culturemods/blob/master/examples/o2_at_1250um.png)


