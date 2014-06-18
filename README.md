A ray tracer made for the final project of my 410 Data Visualization taught by Professor Hank Childs at the University of Oregon. 


It ray traces .VTK files and outputs a png file. Install VTK with CMAKE, and then use CMAKE to make the executable. The files I ray traced were too large to upload to Github. They are currently located at http://www.cs.uoregon.edu/Classes/13F/cis410visualization/final_proj/astro64.vtk and http://www.cs.uoregon.edu/Classes/13F/cis410visualization/final_proj/astro512.vtk. Please note that I do not host those links and don't guarantee they will always be valid. 


Here are some of the images produced with the large datasets of an exploding star, provided by Professor Childs at the links above. 

http://imgur.com/DD3IbvD,Xlcayt8

The first picture is a lower resolution dataset (Astro64.VTK), and the second is the finer resolution (Astro512.VTK)

The code didn't come out perfect but it still produces a rather stunning image of a star (albeit upside down). I believe the issue is in the code that handles matrices.   
