#required downloads: cv2 and send2trash

#!/usr/bin/python3

from picamera2.encoders import H264Encoder, Quality
from picamera2 import Picamera2, MappedArray
import time
import datetime as dt
import os
import cv2
from send2trash import send2trash
from libcamera import controls
from picamera2.outputs import FfmpegOutput

picam2 = Picamera2()

picam2.configure(picam2.create_video_configuration())
picam2.set_controls({'FrameRate': 30}) #specify framerate

#define filename as datestamp
file_name = dt.datetime.now().strftime('%m-%d-%Y_pi#.mp4') #name of video here
file_path = '/home/pi/'
vidname = os.path.join(file_path, file_name)

#define timestamp within the video
colour = (0, 255, 0)
origin = (0, 30)
font = cv2.FONT_HERSHEY_SIMPLEX
scale = 1
thickness = 2
def apply_timestamp(request):
    timestamp = time.strftime("%Y-%m-%d %X")
    with MappedArray(request, "main") as m:
   	 cv2.putText(m.array, timestamp, origin, font, scale, colour, thickness)
picam2.pre_callback = apply_timestamp

#enable encoder
encoder = H264Encoder(10000000)
output = FfmpegOutput(vidname)

picam2.start_recording(encoder, output)
print('recording')
time.sleep(1800) #change this to the length of your video
picam2.stop_recording()

print('uploading')
os.system('rclone sync /home/pi/*.mp4 gdmedia:filename')
send2trash(vidname)
