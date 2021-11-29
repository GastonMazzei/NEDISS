rm FrontUtils/youtube-video.mp4
rm output0.mp4
rm output1.mp4
rm output2.mp4
rm output3.mp4
rm output4.mp4

ffmpeg -loop 1 -i FrontUtils/periodicBackground.png -i graphic/video/result.avi -filter_complex "[1:v]scale=650:-1[fg];[0:v][fg]overlay=(W-w)/2:(H-h)/2:shortest=1" output0.mp4


ffmpeg -i output0.mp4 -filter_complex  \
"[0:v][0:v]overlay=640:0[bg]; \
 [bg][0:v]overlay=640-W,format=yuv420p[out]" \
-map "[out]" -map 0:a:0? -codec:v libx264 -crf 23 -preset medium -c:a copy output1.mp4

ffmpeg -i output1.mp4 -i graphic/video/result-oscillating.avi -filter_complex "[1:v]scale=620:-1[fg];[0:v][fg]overlay=(W-w)/2:(H-h)/2:shortest=1" output2.mp4


ffmpeg -i output2.mp4 -filter_complex  \
"[0:v][0:v]overlay=650:0[bg]; \
 [bg][0:v]overlay=650-W,format=yuv420p[out]" \
-map "[out]" -map 0:a:0? -codec:v libx264 -crf 23 -preset medium -c:a copy output3.mp4

ffmpeg -i output3.mp4  -loop 1  -i graphic/image/maxdistovertime.png -filter_complex "[1:v]scale=600:-1[fg];[0:v][fg]overlay=(W-w)/2:(H-h)/2:shortest=1" output4.mp4

ffmpeg -i output4.mp4 -filter_complex  \
"[0:v][0:v]overlay=610:0[bg]; \
 [bg][0:v]overlay=610-W,format=yuv420p[out]" \
-map "[out]" -map 0:a:0? -codec:v libx264 -crf 23 -preset medium -c:a copy FrontUtils/youtube-video.mp4

rm output0.mp4
rm output1.mp4
rm output2.mp4
rm output3.mp4
rm output4.mp4

