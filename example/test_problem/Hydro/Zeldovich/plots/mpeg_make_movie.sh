ffmpeg -framerate 5/1 -i Combine_Data_%06d.png -c:v libx264 -preset slow -tune animation \
       -pix_fmt yuv420p -s 1502x852 Zeldovich_Collapse_Analysis.mp4