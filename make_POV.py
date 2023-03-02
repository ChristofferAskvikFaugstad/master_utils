from ase.io import write
from ase.io.pov import get_bondpairs
import subprocess
import os
import matplotlib.image as Image
import matplotlib.pyplot as plt
from typing import List
from ase.atoms import Atoms

povray_path = "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"

radii_dict = {
        'Ni': 0.8,
        'H': 0.2,
        "C" : 0.5,
        "O" : 0.6,
    }

class POVMaker:
    def __init__(self, images : List[Atoms], folder : str):
        self.images = images
        self.folder = folder
        os.makedirs(folder, exist_ok=True)
    

    def make_png(self, i,
        radii_dict = radii_dict,
        bond_radius = 1.1,
        show_unit_cell = 0,
        rotation = '33x,15y,10z',
        camera_dist = 25,
        canvas_width = 400,
        move_center = (0,0,0),
        **kwargs
        ):
        image = self.images[i]
        image.translate(move_center)
        r = [radii_dict[atom.symbol] for atom in image]
        bondpairs = get_bondpairs(image, radius=bond_radius)
        folder = self.folder
        filename = f'{folder}/{i}.pov'
        config = f'{i}.ini'
        write(filename, image, format='pov',
                        radii=r, rotation=rotation,show_unit_cell = show_unit_cell,
                        povray_settings=dict(canvas_width=canvas_width, bondatoms=bondpairs, camera_dist = camera_dist, **kwargs))
        
        cwd = os.getcwd()
        os.chdir(folder)
        subprocess.run(["pvengine64.exe", "/EXIT", "/RENDER", f"{config}"])
        os.remove(f"{i}.ini")
        os.remove(f"{i}.pov")
        os.chdir(cwd)
    
    def make_pngs(self,
        radii_dict = radii_dict,
        bond_radius = 1.1,
        show_unit_cell = 0,
        rotation = '33x,15y,10z',
        camera_dist = 25,
        canvas_width = 400,
        **kwargs
        ):
        for i,image in enumerate(self.images):
            self.make_png(i,radii_dict,bond_radius,show_unit_cell, rotation,camera_dist,canvas_width,**kwargs)
    
    
    def preview(self, i, 
        radii_dict = radii_dict,
        bond_radius = 1.1,
        show_unit_cell = 0,
        rotation = '33x,15y,10z',
        camera_dist = 25,
        canvas_width = 400,
        **kwargs):
        
        self.make_png(i,radii_dict,bond_radius,show_unit_cell, rotation,camera_dist,canvas_width,**kwargs)
        
        png = Image.imread(f"{self.folder}/{i}.png")
        plt.imshow(png)

    def make_gif(self, 
        save_all : bool = True,
        duration : int = 300,
        loop : int = 0,
         **kwargs):
            """
            In the cwd it gatters all png files and sorts on the float of the
            base name of the pngs to order
            """    
            from PIL import Image
            import glob

            frames =  []
            imgs = glob.glob(f"{self.folder}/*.png")

            imgs = sorted(imgs, key = lambda x: float(os.path.basename(x)[:-4]) )

            for i in imgs:
                new_frame = Image.open(i)
                frames.append(new_frame)

            frames[0].save(f"{self.folder}/gif.gif", format='GIF',
                            append_images = frames[1:],
                            save_all = save_all,
                            duration = duration, loop = loop, **kwargs)





if __name__ == "__main__":
    from my_modules.project import *
    path = "H2_adsorbtion_neb/neb_os_s/neb.traj"
    traj = Trajectory(path)

    images = [image for image in traj]
    
    pov = POVMaker(images, "test_POV")

    # pov.make_png(1)

    pov.preview(2, move_center = (10,0,0))

    plt.show()