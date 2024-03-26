import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read the data
GC_path = '' 
Center_path = ''
Start_Time = 2000
End_Time = 5000
x_limit = (7000,9000)
y_limit = (7000,9000)

Center_df = pd.read_csv(Center_path, sep='\s+', header=0)
GC_df = pd.read_csv(GC_path, sep='\s+', header=0)
Center_df = Center_df[(Center_df['Time'] > Start_Time) & (Center_df['Time'] < End_Time)]
GC_df = GC_df[(GC_df['Time'] > Start_Time) & (GC_df['Time'] < End_Time)]

# plot the trajectory of GC and center with different colors at different time and make a movie
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlabel('X (pc)')
ax.set_ylabel('Y (pc)')
ax.set_xlim(x_limit[0], x_limit[1])
ax.set_ylim(y_limit[0], y_limit[1])
ax.scatter(Center_df['CenterX'],Center_df['CenterY'],color='b',label='Center')
line_center, = ax.plot([], [], color='b', label='Center', linewidth=0.5)
line_gc, = ax.plot([], [], color='r', label='GC', linewidth=0.75)
ax.legend(loc='upper right')

def animate(i):
    line_center.set_data(Center_df['CenterX'][i-500:i], Center_df['CenterY'][i-500:i])
    line_gc.set_data(GC_df['GCX'][i-500:i], GC_df['GCY'][i-500:i])
    current_time = Center_df['Time'].iloc[i]
    ax.set_title(f'Trajectory of GC and Center at Time {current_time}')
    return line_center,line_gc

# Determine the number of frames for the animation
num_frames = min(len(Center_df), len(GC_df))

print("Making animation...")
ani = FuncAnimation(fig, animate, frames=num_frames, interval=20, blit=True)

# Save the animation
ani.save('gc_center_trajectory.mp4', writer='ffmpeg')
print("Making animation...Done")


