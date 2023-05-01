#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : batch_rename.py
@Time     : 2023/05/02 00:05:54
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib

import os


def batch_rename(file_dir, suffix):
    """
    Batch modify the names of all directories name under the specified directory.

    Args:
    - file_dir (str): The path of the directory to be modified.
    - suffix (str): The suffix to be added to or removed from the directory names.

    Returns:
    - None
    """
    for root, dirs, files in os.walk(file_dir):
        for dir in dirs:
            try:
                new_dir = dir if suffix in dir else dir + suffix
                os.rename(os.path.join(root, dir), os.path.join(root, new_dir))
            except Exception as e:
                print(f'\nWhen renaming the folder: {dir}, an error occurred: {e}!')
                continue


def file_rename(file_dir, suffix):
    """
    Batch modify the names of all files under the specified directory.

    Args:
    - file_dir (str): The path of the directory to be modified.
    - suffix (str): The suffix to be added to or removed from the file names.

    Returns:
    - None
    """
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if file.endswith(".tsv"):
                try:
                    folder_name, _ = os.path.splitext(
                        os.path.basename(root)
                    )  # Exp112121_F1_R2
                    folder_name = folder_name.replace(suffix, '')  # Exp112121
                    new_file = file.replace(
                        folder_name, ''
                    )  # Exp112121_Colon_049761.tsv
                    while suffix in new_file or suffix in folder_name:
                        new_file = new_file.replace(suffix, '')
                        folder_name = folder_name.replace(suffix, '')
                        new_file = new_file.replace(folder_name, '')
                    new_file = f"{folder_name}{suffix}_{new_file}"

                    if new_file != file:
                        os.rename(
                            os.path.join(root, file), os.path.join(root, new_file)
                        )
                except Exception as e:
                    print(f'\nWhen renaming the file: {file}, an error occurred: {e}!')
                    continue


def batch_renaming(file_dir, suffix):
    file_rename(file_dir, suffix)
    batch_rename(file_dir, suffix)


if __name__ == "__main__":
    file_dir = r"/Users/ranpeng/Desktop/DIA-QC/file_transfer/sichuan/"
    suffix = "_F1_R2"
    batch_renaming(file_dir, suffix)
