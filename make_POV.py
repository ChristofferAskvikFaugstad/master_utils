from ase.io import write
from ase.io.pov import get_bondpairs
import subprocess
import os
import matplotlib.image as Image
import matplotlib.pyplot as plt
from typing import List
from ase.atoms import Atoms
import numpy as np
import matplotlib as mpl


povray_path = "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"

radii_dict = {
    "Ni": 0.8,
    "H": 0.2,
    "C": 0.5,
    "O": 0.6,
}


class POVMaker:
    @property
    def get_default_cmap(self):
        return mpl.cm.get_cmap("RdBu")

    def __init__(self, images: List[Atoms], folder: str):
        self.images = images
        self.folder = folder
        os.makedirs(folder, exist_ok=True)

    def make_png(
        self,
        i,
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        move_center=(0, 0, 0),
        **kwargs,
    ):
        image = self.images[i]
        image.translate(move_center)
        r = [radii_dict[atom.symbol] for atom in image]
        bondpairs = get_bondpairs(image, radius=bond_radius)
        folder = self.folder
        filename = f"{folder}/{i}.pov"
        config = f"{i}.ini"
        write(
            filename,
            image,
            format="pov",
            radii=r,
            rotation=rotation,
            show_unit_cell=show_unit_cell,
            povray_settings=dict(
                canvas_width=canvas_width,
                bondatoms=bondpairs,
                camera_dist=camera_dist,
                **kwargs,
            ),
        )

        cwd = os.getcwd()
        os.chdir(folder)
        subprocess.run(["pvengine64.exe", "/EXIT", "/RENDER", f"{config}"])
        os.remove(f"{i}.ini")
        os.remove(f"{i}.pov")
        os.chdir(cwd)

    def make_pngs(
        self,
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        **kwargs,
    ):
        for i, image in enumerate(self.images):
            self.make_png(
                i,
                radii_dict,
                bond_radius,
                show_unit_cell,
                rotation,
                camera_dist,
                canvas_width,
                **kwargs,
            )

    def preview(
        self,
        i,
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        **kwargs,
    ):

        self.make_png(
            i,
            radii_dict,
            bond_radius,
            show_unit_cell,
            rotation,
            camera_dist,
            canvas_width,
            **kwargs,
        )

        png = Image.imread(f"{self.folder}/{i}.png")
        plt.imshow(png)

    def make_gif(
        self,
        save_all: bool = True,
        duration: int = 300,
        loop: int = 0,
        name: str = "gif",
        image_name="",
        **kwargs,
    ):
        """
        In the cwd it gatters all png files and sorts on the float of the
        base name of the pngs to order
        """
        from PIL import Image
        import glob

        frames = []
        imgs = glob.glob(f"{self.folder}/*.png")

        imgs = sorted(imgs, key=lambda x: os.path.basename(x)[:-4])

        for i in imgs:
            new_frame = Image.open(i)
            frames.append(new_frame)

        frames[0].save(
            f"{self.folder}/{name}{image_name}.gif",
            format="GIF",
            append_images=frames[1:],
            save_all=save_all,
            duration=duration,
            loop=loop,
            **kwargs,
        )

    def make_induvidual_color(
        self,
        i,
        colors,
        name="C",
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        move_center=(0, 0, 0),
        **kwargs,
    ):
        image = self.images[i]
        image.translate(move_center)
        r = [radii_dict[atom.symbol] for atom in image]
        bondpairs = get_bondpairs(image, radius=bond_radius)
        folder = self.folder
        filename = f"{i}{name}.pov"
        config = f"{i}{name}.ini"

        cwd = os.getcwd()
        os.chdir(folder)
        write(
            filename,
            image,
            format="pov",
            radii=r,
            rotation=rotation,
            show_unit_cell=show_unit_cell,
            colors=colors,
            povray_settings=dict(
                canvas_width=canvas_width,
                bondatoms=bondpairs,
                camera_dist=camera_dist,
                **kwargs,
            ),
        )

        subprocess.run(["pvengine64.exe", "/EXIT", "/RENDER", f"{config}"])
        os.remove(config)
        os.remove(filename)
        os.chdir(cwd)

    def get_colors(self, array, norm, cmap=None):
        """
        Returns a list of colors from a list of values
        """
        if cmap is None:
            cmap = self.get_default_cmap
        return [cmap(norm(i)) for i in array]

    def get_attribute_array(self, func):
        """
        Makes an array of the attribute of the atoms in the images
        """
        array = []
        for image in self.images:
            array += list(func(image))
        # array = np.array([func(image) for image in self.images]).flatten()
        array = np.array(array).flatten()
        return array

    def get_norm(self, array, center = False):
        """
        Returns a normed array
        """
        if center:
            return mpl.colors.TwoSlopeNorm(vcenter = 0, vmin = np.min(array),vmax = np.max(array))
        else:
            return mpl.colors.Normalize(np.min(array), np.max(array))

    def make_colorbar(self, array, cmap=None, orientation="horizontal", label="", nticks=5, center = False):
        """
        Makes a colorbar from a list of colors
        """
        if cmap is None:
            cmap = self.get_default_cmap

        fig = plt.figure()
        if orientation == "horizontal":
            ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])
        if orientation == "vertical":
            ax = fig.add_axes([0.80, 0.05, 0.05, 0.45])
        if center:
            cb = mpl.colorbar.ColorbarBase(
            ax,
            orientation=orientation,
            cmap=cmap,
            norm=self.get_norm(array, center = center),  # vmax and vmin
            ticks=[*np.linspace(np.min(array), 0, int(nticks/2), endpoint=False),*np.linspace(0,np.max(array), int(nticks/2)+1, endpoint=True)],
            label=label,
            )
        
        else:
            cb = mpl.colorbar.ColorbarBase(
            ax,
            orientation=orientation,
            cmap=cmap,
            norm=self.get_norm(array),  # vmax and vmin
            ticks=np.linspace(np.min(array), np.max(array), nticks),
            label=label,
        )

        plt.savefig(f"{self.folder}/just_colorbar", bbox_inches="tight")

    def make_png_attributes(
        self,
        i,
        func,
        name="C",
        cmap=None,
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        center = False,
        **kwargs,
    ):
        array = self.get_attribute_array(func)
        norm = self.get_norm(array, center = center)
        image = self.images[i]

        cur_array = np.array([func(image)]).flatten()
        colors = self.get_colors(cur_array, norm, cmap=cmap)
        self.make_induvidual_color(
            i,
            colors,
            name=name,
            radii_dict=radii_dict,
            bond_radius=bond_radius,
            show_unit_cell=show_unit_cell,
            rotation=rotation,
            camera_dist=camera_dist,
            canvas_width=canvas_width,
            **kwargs,
        )

    def make_pngs_attributes(
        self,
        func,
        name="C",
        cmap=None,
        radii_dict=radii_dict,
        bond_radius=1.1,
        show_unit_cell=0,
        rotation="33x,15y,10z",
        camera_dist=25,
        canvas_width=400,
        center = False,
        **kwargs,
    ):
        array = self.get_attribute_array(func)
        norm = self.get_norm(array, center = center)

        for i, image in enumerate(self.images):
            cur_array = np.array([func(image)]).flatten()
            colors = self.get_colors(cur_array, norm, cmap=cmap)
            self.make_induvidual_color(
                i,
                colors,
                name=name,
                radii_dict=radii_dict,
                bond_radius=bond_radius,
                show_unit_cell=show_unit_cell,
                rotation=rotation,
                camera_dist=camera_dist,
                canvas_width=canvas_width,
                **kwargs,
            )

    def make_horisontal_figure(
        self,
        save_all: bool = True,
        duration: int = 300,
        loop: int = 0,
        name: str = "gif",
        image_name="",
        **kwargs,
    ):
        """
        In the cwd it gatters all png files and sorts on the float of the
        base name of the pngs to order
        """
        from PIL import Image
        import glob

        frames = []
        imgs = glob.glob(f"{self.folder}/*.png")

        imgs = sorted(imgs, key=lambda x: os.path.basename(x)[:-4])

        for i in imgs:
            new_frame = Image.open(i)
            frames.append(new_frame)

        frames[0].save(
            f"{self.folder}/{name}{image_name}.gif",
            format="GIF",
            append_images=frames[1:],
            save_all=save_all,
            duration=duration,
            loop=loop,
            **kwargs,
        )



def isosurface():
    from ase.calculators.vasp import VaspChargeDensity
    from ase.io import write

    spin_cut_off = 0.4
    density_cut_off = 0.00043

    # rotation = '24x, 34y, 14z'
    rotation = '90x,90y,0z'
    rotation = '0x,0y,0z'

    name = "upLUMO"
    vchg = VaspChargeDensity("in/song_es_full_test2/upLUMO.vasp")
    atoms = vchg.atoms[0]

    povray_settings = {
        # For povray files only
        'pause': False,  # Pause when done rendering (only if display)
        'transparent': False,  # Transparent background
        'canvas_width': 400,  # Width of canvas in pixels
        'canvas_height': None,  # Height of canvas in pixels
        'camera_dist': 25.0,  # Distance from camera to front atom
        # 'camera_type': 'orthographic angle 35',  # 'perspective angle 20'
        'textures': len(atoms) * ['ase3']}

    # some more options:
    # 'image_plane'  : None,  # Distance from front atom to image plane
    #                         # (focal depth for perspective)
    # 'camera_type'  : 'perspective', # perspective, ultra_wide_angle
    # 'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
    # 'area_light'   : [(2., 3., 40.) ,# location
    #                   'White',       # color
    #                   .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
    # 'background'   : 'White',        # color
    # 'textures'     : tex, # Length of atoms list of texture names
    # 'celllinewidth': 0.05, # Radius of the cylinders representing the cell

    # atoms.edit()
    generic_projection_settings = {
        'rotation': rotation,
        'radii': atoms.positions.shape[0] * [0.3],
        # 'show_unit_cell': 1
        }

    # write returns a renderer object which needs to have the render method called

    write(f'{name}.pov', atoms,
        **generic_projection_settings,
        povray_settings=povray_settings,
        isosurface_data=dict(density_grid=vchg.chg[0],
                            cut_off=density_cut_off,
                            color=[0.25, 0.25, 0.80, 0.1],))

    subprocess.run(["pvengine64.exe", "/EXIT", "/RENDER", f"{name}.ini"])
    # spin up density, how to specify color and transparency r,g,b,t and a
    # material style from the standard ASE set


    # write('NiO_marching_cubes2.pov', atoms,
    #     **generic_projection_settings,
    #     povray_settings=povray_settings,
    #     isosurface_data=dict(density_grid=vchg.chgdiff[0],
    #                         cut_off=density_cut_off,
    #                         closed_edges=True,
    #                         color=[0.25, 0.25, 0.80, 0.1],
    #                         material='simple')).render()

    # # spin down density, how to specify a povray material
    # # that looks like pink jelly
    # fun_material = '''
    # material {
    #     texture {
    #     pigment { rgbt < 0.8, 0.25, 0.25, 0.5> }
    #     finish{ diffuse 0.85 ambient 0.99 brilliance 3 specular 0.5 \
    # roughness 0.001
    #         reflection { 0.05, 0.98 fresnel on exponent 1.5 }
    #         conserve_energy
    #     }
    #     }
    #     interior { ior 1.3 }
    # }
    # photons {
    #     target
    #     refraction on
    #     reflection on
    #     collect on
    # }'''

    # write('NiO_marching_cubes3.pov', atoms,
    #     **generic_projection_settings,
    #     povray_settings=povray_settings,
    #     isosurface_data=dict(density_grid=vchg.chgdiff[0],
    #                         cut_off=-spin_cut_off,
    #                         gradient_ascending=True,
    #                         material=fun_material)).render()


if __name__ == "__main__":
    from utils.master import *

    text = """ 1       -0.005  -0.006   0.807   0.797
        2       -0.006  -0.010   0.766   0.751
        3       -0.003  -0.014   0.721   0.704
        4       -0.003  -0.014   0.720   0.703
        5        0.016  -0.017   0.763   0.762
        6       -0.003  -0.009   0.684   0.671
        7       -0.004  -0.010   0.685   0.672
        8       -0.012  -0.018   0.728   0.698
        9       -0.012  -0.018   0.729   0.699
    10       -0.007  -0.015   0.654   0.631
    11       -0.001  -0.012   0.743   0.730
    12       -0.001  -0.012   0.743   0.731
    13       -0.013  -0.012   0.452   0.426
    14       -0.013  -0.013   0.454   0.429
    15       -0.009  -0.014   0.764   0.741
    16       -0.009  -0.014   0.765   0.741
    17       -0.009  -0.023   0.751   0.719
    18       -0.009  -0.023   0.756   0.724
    19       -0.011  -0.027   0.763   0.726
    20       -0.011  -0.027   0.763   0.726
    21       -0.012  -0.027   0.778   0.739
    22       -0.012  -0.027   0.778   0.739
    23       -0.014  -0.023   0.425   0.388
    24       -0.007  -0.026   0.750   0.718
    25       -0.007  -0.026   0.749   0.717
    26       -0.007  -0.029   0.804   0.767
    27       -0.009  -0.025   0.765   0.731
    28       -0.007  -0.024   0.679   0.649
    29       -0.007  -0.024   0.679   0.649
    30       -0.007  -0.022   0.568   0.539
    31       -0.001  -0.050   0.000  -0.051
    32        0.002  -0.022   0.000  -0.020"""
    text = """    1        0.499   0.212   8.376   9.086
    2        0.484   0.262   8.376   9.122
    3        0.478   0.263   8.385   9.127
    4        0.479   0.264   8.385   9.128
    5        0.500   0.267   8.365   9.132
    6        0.494   0.272   8.385   9.152
    7        0.494   0.273   8.385   9.153
    8        0.491   0.305   8.373   9.169
    9        0.491   0.305   8.372   9.168
   10        0.476   0.254   8.389   9.119
   11        0.495   0.289   8.372   9.156
   12        0.495   0.290   8.371   9.155
   13        0.460   0.388   8.426   9.274
   14        0.460   0.388   8.425   9.273
   15        0.475   0.320   8.346   9.141
   16        0.476   0.321   8.346   9.143
   17        0.492   0.337   8.354   9.183
   18        0.491   0.339   8.352   9.182
   19        0.495   0.381   8.334   9.211
   20        0.496   0.383   8.334   9.213
   21        0.493   0.389   8.335   9.218
   22        0.493   0.390   8.335   9.219
   23        0.497   0.500   8.407   9.404
   24        0.489   0.380   8.339   9.208
   25        0.489   0.380   8.340   9.208
   26        0.500   0.392   8.327   9.219
   27        0.543   0.605   8.349   9.498
   28        0.519   0.569   8.342   9.430
   29        0.519   0.568   8.342   9.429
   30        0.509   0.561   8.346   9.417
   31        1.631   3.503   0.000   5.134
   32        0.917   1.597   0.000   2.514"""

    isosurface()

    # magnetisaiton = np.array([float(line.split()[-1]) for line in text.split("\n")])
    # path = "in/song_es"
    # structure = get_structure_path(path, relative = False)
    # bader_charges = get_bader_charges(path, relative = False, valency=[10]*30)
    # structure.info["bader_charges"] = bader_charges

    # def func(image):
    #     return image.info["bader_charges"]
    
    # structures = [structure.copy() for _ in range(3)]
    # structures[1].rotate(90,"x")
    # structures[2].rotate(90,"y")
    # pov = POVMaker(structures, "test_POV")
    # pov.make_pngs_attributes(func, rotation="0x,0y,0z", center=True)
    # array = pov.get_attribute_array(func)
    # pov.make_colorbar(array, orientation="vertical", label="Bader charge [e]", center=True, nticks = 6)


    # structure = get_structure_path("single/ni30_COs/CO12-13-22")
    # norm = mpl.colors.Normalize(np.min(magnetisaiton), np.max(magnetisaiton))
    # print(norm(magnetisaiton))

    # structures = get_all_structures_path("single/ni30_COs/CO12-13-22")[:2]
    # pov.make_pngs_attributes(
    #     Atoms.get_magnetic_moments,
    #     rotation="180x,0y,180z",
    #     cmap=mpl.cm.get_cmap("Reds"),
    # )
    # for folder in get_vasp_calculation_names_path(
    #     "neb/CO_dis/CO5-25_C2-5-12-25O6-13-25"
    # ):
    #     print(folder)
    #     structure = get_structure_path(
    #         os.path.join("neb/CO_dis/CO5-25_C2-5-12-25O6-13-25", folder)
    # )



    # images = get_structures_folder("single/ni30_2COs")
    # pov = POVMaker(images, "test_POV")
    # pov.make_colorbar(magnetisaiton, orientation="vertical")

    # pov.make_pngs_attributes(Atoms.get_magnetic_moments, rotation="46x,-77y,-39z")
    # pov.make_gif(name="magmoms", image_name="C", rotation="46x,-77y,-39z")

    # images = [structure]
    # pov = POVMaker(images, "test_POV2")
    # pov.make_colorbar(magnetisaiton)

    # import matplotlib

    # cmap = matplotlib.cm.get_cmap("Blues")

    # colors = cmap(
    #     (magnetisaiton - np.min(magnetisaiton))
    #     / (np.max(magnetisaiton) - np.min(magnetisaiton))
    # )
    # import matplotlib.pyplot as plt
    # import matplotlib as mpl

    # fig = plt.figure()
    # ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])

    # cb = mpl.colorbar.ColorbarBase(
    #     ax,
    #     orientation="horizontal",
    #     cmap="Blues",
    #     norm=mpl.colors.Normalize(
    #         np.min(magnetisaiton), np.max(magnetisaiton)
    #     ),  # vmax and vmin
    #     label="Magnetisation",
    #     ticks=np.linspace(np.min(magnetisaiton), np.max(magnetisaiton), 5),
    # )

    # plt.savefig("test_POV2/just_colorbar", bbox_inches="tight")
    # # view(structure, block=True)

    # import numpy as np

    # pov.make_induvidual_color(0, colors, rotation="180x,0y,180z")

    # from PIL import Image

    # list_im = ["test_POV2/0.png", "test_POV2/just_colorbar.png"]
    # imgs = [Image.open(i) for i in list_im]

    # png = imgs[0]
    # colorbar = imgs[1]

    # colorbar = colorbar.resize(
    #     (png.size[0], int(png.size[0] / colorbar.size[0] * colorbar.size[1]))
    # )

    # imgs_comb = np.vstack([png, colorbar])
    # imgs_comb = Image.fromarray(imgs_comb)
    # imgs_comb.save("test_POV2/0_colorbar.png")
