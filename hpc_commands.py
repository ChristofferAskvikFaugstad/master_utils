import shutil
import os


def copy_neb_files(destination: str, num_images: int):
    """
    Copies the files needed for a NEB calculation to a new directory.
    """
    # Make the new directory
    os.makedirs(destination, exist_ok=True)
    for i in range(num_images + 2):
        os.makedirs(os.path.join(destination, f"{i:02}"), exist_ok=True)
    # Copy the files
    for i in [0, num_images + 1]:
        shutil.copyfile(
            os.path.join(os.getcwd(), f"{i:02}", "vasprun.xml"),
            os.path.join(destination, f"{i:02}", "vasprun.xml"),
        )
    for i in range(1, num_images + 1):
        shutil.copyfile(
            os.path.join(os.getcwd(), f"{i:02}", "vasprun.xml"),
            os.path.join(destination, f"{i:02}", "vasprun_prev.xml"),
        )
        shutil.copyfile(
            os.path.join(os.getcwd(), f"{i:02}", "OUTCAR"),
            os.path.join(destination, f"{i:02}", "OUTCAR_prev.xml"),
        )
        shutil.copyfile(
            os.path.join(os.getcwd(), f"{i:02}", "CONTCAR"),
            os.path.join(destination, f"{i:02}", "POSCAR"),
        )


copy_neb_files("/cluster/home/chrisafa/HCOtoCH_O_neb/job", 5)


text = """slurm-8484639_0.out:Starting job 8487962 on c1-3,c3-2,c5-14,c10-47 on saga at Thu Apr 6 03:30:38 CEST 2023        
slurm-8484639_1.out:Starting job 8488682 on c1-41,c2-44,c3-9,c5-28 on saga at Thu Apr 6 04:45:42 CEST 2023        
slurm-8484639_2.out:Starting job 8488818 on c1-15,c10-[27-28],c11-16 on saga at Thu Apr 6 05:34:34 CEST 2023      
slurm-8484639_3.out:Starting job 8488830 on c1-41,c2-44,c3-9,c5-28 on saga at Thu Apr 6 07:07:01 CEST 2023        
slurm-8484639_4.out:Starting job 8488840 on c1-3,c3-2,c5-14,c10-47 on saga at Thu Apr 6 07:51:45 CEST 2023        
slurm-8484639_5.out:Starting job 8484639 on c1-15,c10-[27-28],c11-16 on saga at Thu Apr 6 08:54:57 CEST 2023"""

ints = []
ids = []
for line in text.split("\n"):
    ints.append(int(float(line.split("_")[1][0:2])))
    ids.append(int(line.split(" ")[2]))
print(ints)
print(ids)


import shutil
import os
def copy_slurm_ID_task_output(task_IDs, job_IDs, output_dir):
    for task_ID, job_ID in zip(task_IDs, job_IDs):
        calculaiton_dir = os.path.join("/cluster/work/users/chrisafa/temp/", str(job_ID), str(task_ID))
        cur_output_dir = os.path.join(output_dir, str(task_ID))
        os.makedirs(cur_output_dir, exist_ok=True)
        for file in ["OUTCAR","vasprun.xml", "OSZICAR"]:
            shutil.copyfile(os.path.join(calculaiton_dir, file), os.path.join(cur_output_dir, file))

task_IDs = [0, 1, 2, 3, 4, 5]
job_IDs = [8487962, 8488682, 8488818, 8488830, 8488840, 8484639]
output_dir = "/cluster/home/chrisafa/out/CH_O_scr"
copy_slurm_ID_task_output(task_IDs, job_IDs, output_dir)
