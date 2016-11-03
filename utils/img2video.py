# Copyright (C) 2016 Zhixian MA <zxma_sjtu@qq.com>
# MIT license

"""
A simple tool to translate a group of images to a video. In the gas track work,
it is used to gather the simulated gas images into a portable video.

References
----------
[1] OpenCV-Python Tutorials
    https://opencv-python-toturials.readthedocs.io/en/latest/
[2] Solem, J.E.
    "Programming Computer Vision with Python"
    O'Reilly, 2012.
"""

import cv2
import os
import re

def img2video(input_dir, fileout, fps=4.0,imgwidth=(800,800),
              fxpression = "snap"):
    """
    Image to video transformation.

    Parameters
    ----------
    input_dir: string
        Name of the folder holding those images.
    fxpression: string
        The regular expression of those images.
    fileout: string
        The filepath of output.
    fps: float
        Frames per second.
    imgwidth: tuple
        Width and height of the video frame.

    Reference
    ---------
    cv2.VideoWriter
    https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/
    py_gui/py_video_display/py_video_display.html?highlight=videowriter
    """

    # Define the codec and create ViderWriter object
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    output = cv2.VideoWriter(fileout,fourcc,fps,imgwidth)

    # Detect pngfiles
    files = os.listdir(input_dir)
    files.sort()
    for f in files:
        if len(re.findall(r"snap",f)):
            print(f)
            image = cv2.imread(f)
            output.write(image)

    # Release
    output.release()
