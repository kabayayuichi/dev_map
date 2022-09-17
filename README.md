# dev_map
A python program drawing the developing map of a surface group representation to PSL(2,C) parametrized by Fenchel-Nielsen coordinates.

## How to use?
Download dev_map.py and demo_g2.py. Then type "python demo_g2.py". It creats two eps files "demo_g2.eps", "demo_g2_grafting.eps", and two povray files "demo_g2.pov", "demo_g2_ball.pov".
The eps files are too large, so it might be better to convert into jpeg files using e.g. ImageMagik ("convert demo_g2.eps demo_g2.jpg"). 
To convert povray files, use e.g. "povray -d demo_g2.pov".

You can adjast Fenchel-Nielsen parameters in demo_g2.py. But if the surface group representation is indiscete, eps and pov files will be too large. I recommend set the number "num_iteration" in demo_g2.py to 5~8 on the first try.
