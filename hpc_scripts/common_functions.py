import os
def find_numbered_folders(directory):
    numbered_folders = []

    for folder_name in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, folder_name)):
            if folder_name.isdigit():
                numbered_folders.append(folder_name)
    numbered_folders.sort()
    return numbered_folders

