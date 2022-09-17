# dev_map
A python program drawing the developing map of a surface group representation to PSL(2,C) parametrized by Fenchel-Nielsen coordinates.

## How to use?
Download dev_map.py and demo_g2.py. Then type "python demo_g2.py". It creates two eps files "demo_g2.eps", "demo_g2_grafting.eps", and two povray files "demo_g2.pov", "demo_g2_ball.pov".
The eps files are too large (40~50MB), so it might be better to convert into jpeg files using e.g. ImageMagik ("convert demo_g2.eps demo_g2.jpg"). 
To convert povray files, use e.g. "povray -d demo_g2.pov".

You can adjust Fenchel-Nielsen parameters in demo_g2.py. But if the surface group representation is indiscete, eps and pov files will be too much large (maybe over 10GB). I recommend set the number "num_iteration" in demo_g2.py to 5~8 on the first try.

<img src="https://user-images.githubusercontent.com/38853563/190839133-251e1628-c3cf-4aba-9157-6da953625e95.jpg" width="300px"><img src="https://user-images.githubusercontent.com/38853563/190839139-73ee89ab-5039-400f-b9c3-3fad3ce9a02c.jpg" width="300px">
<img src="https://user-images.githubusercontent.com/38853563/190839141-3a145403-e1d0-4529-8a3a-6bd0844fe712.png" width="300px">
<img src="https://user-images.githubusercontent.com/38853563/190839144-b53962f1-930e-4f91-823a-b236f0699487.png" width="300px">
