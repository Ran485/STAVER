#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : download_pride_data.py
@Time     : 2022/12/06 21:15:33
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2019-2020, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib

import ppx
import pandas as pd
import os
import requests
from multiprocessing.pool import ThreadPool
from rich.progress import track
import datetime
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import cpu_count

# from utils import create_output_dir

# Get the number of CPU threads
NUM_WORKERS = cpu_count() - 2


def create_output_dir(out_path=None, creat_time_subdir=False):
    """create the output directory

    Args:
    -------------
        filename ([str]): A given filename

        creat_time_subdir (bool, optional):
                        creat 2021-11-12 subdirectory,defaults to False.
    Returns:
    -------------
        output_dir [str]: the output directory
    """
    if not os.path.exists(out_path):  # in case root_dir doesn't exist
        os.makedirs(out_path)
        print("Successfully created output subdir: {}".format(out_path))
        if creat_time_subdir:
            date_string = datetime.now().strftime("%Y_%m_%d")
            out_path = os.path.join(out_path, date_string)
            if not os.path.isdir(out_path):
                os.mkdir(out_path)
                print("Successfully created output subdir: {}".format(out_path))
    else:
        print("The current path: {} already exist".format(out_path))
    return out_path + "/"


def find_project_files(Pride_ID, outpath) -> list:
    global proj
    proj = ppx.find_project(Pride_ID, local=outpath)
    # Print remote files before downloading
    print("Files to be downloaded:", proj.remote_files())
    return proj, proj.remote_files()


def find_small_file(filepath):
    try:
        size = os.path.getsize(filepath)
        file_size = 2 * 1024
        if size > file_size:
            print(f"\nThe file {filepath} size is : {size}")
            return None
        else:
            return filepath
    except Exception as e:
        print(f"e")


def conut_empty_file(filepath):
    try:
        empty_nums = 0
        for maindir, subdir, file_name_list in os.walk(filepath):
            for filename in file_name_list:
                file = os.path.join(maindir, filename)
                size = os.path.getsize(file)
                file_size = 2 * 1024
                if size < file_size:
                    empty_nums += 1
    except Exception as e:
        print(f"e")
    return empty_nums


def download_data(filename):
    proj.download(filename)


def validate_file(Pride_ID, outpath):
    _, files = find_project_files(Pride_ID, outpath)
    filepaths = [(outpath + file) for file in files]
    return filepaths


def single_threading(Pride_ID, outpath):
    filepaths = validate_file(Pride_ID, outpath)
    project, filenames = find_project_files(Pride_ID=Pride_ID, outpath=outpath)
    for filename in filenames:
        filepath = outpath + filename
        if not os.path.exists(filepath):
            download_data(filename)


def multi_threading(Pride_ID, outpath, thread=NUM_WORKERS):
    proj, files = find_project_files(Pride_ID=Pride_ID, outpath=outpath)
    filenames = []
    for file in files:
        filepath = outpath + file
        if not os.path.exists(filepath):
            filenames.append(file)
        else:
            filename = find_small_file(filepath)
            if filename != None:
                filenames.append(file)
    with ThreadPoolExecutor(thread) as threads:
        res = threads.map(download_data, filenames)


if __name__ == "__main__":
    # download data from firmiana
    PRIDE_ID = "PXD030304"
    global outpath
    # outpath = create_output_dir('/Volumes/RanPeng_Pas/肿瘤新抗原/Raw_Data/Jiang_Nature')
    outpath = r"/Volumes/Disk3/STAVER-revised/"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    else:
        print(f"The current path: {outpath} already exist!\n")
    print(f"Calling {NUM_WORKERS} threads to download data!\n")
    multi_threading(PRIDE_ID, outpath, thread=NUM_WORKERS)
    empty_nums = conut_empty_file(outpath)
    print(empty_nums)
    while empty_nums > 0:
        multi_threading(PRIDE_ID, outpath, thread=NUM_WORKERS)  ## 多线程并行下载，默认20个线程
        empty_nums = conut_empty_file(outpath)
    # single_threading(PRIDE_ID, outpath) ## 单线程下载

    # proj = ppx.find_project(
    #     "PXD030304", local="/Volumes/Disk3/STAVER-revised/PXD030304"
    # )
    # print(proj.remote_files())
# proj.download("Plasma_1-B10_22_304.wiff.scan")
