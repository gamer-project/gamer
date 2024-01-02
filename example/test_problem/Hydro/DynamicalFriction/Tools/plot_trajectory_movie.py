import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Feed in the data
GC_path = '../Record__RESULT_GC_position.txt' 
Center_path = '../Record__RESULT_Center.txt'

Center_df = pd.read_csv(Center_path, sep='\s+', header=0)
GC_df = pd.read_csv(GC_path, sep='\s+', header=0)
Center_df = Center_df[(Center_df['Time'] > 2000) & (Center_df['Time'] < 5000)]
GC_df = GC_df[(GC_df['Time'] > 2000) & (GC_df['Time'] < 5000)]

print(Center_df.head())

#plt.plot(Center_df['CenterX'], Center_df['CenterY'], color='b', label='Center', linewidth=0.5)
#plt.plot(GC_df['GCX'], GC_df['GCY'], color='r', label='GC', linewidth=0.75)
#plt.xlabel('X (pc)')
#plt.ylabel('Y (pc)')
#plt.title('Trajectory of GC and Center')
#plt.legend(loc='upper right')
#plt.savefig('trajectory.png')
#exit()

# plot the trajectory of GC and center with different colors at different time and make a movie
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlabel('X (pc)')
ax.set_ylabel('Y (pc)')
ax.set_xlim(7000, 9000)
ax.set_ylim(7000, 9000)
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


# Center_df = Center_df[(Center_df['Time'] > 4500) & (Center_df['Time'] < 8000)]

# fig, ax = plt.subplots(figsize=(10, 10))
# ax.set_xlabel('X (pc)')
# ax.set_ylabel('Y (pc)')
# ax.set_xlim(12000, 13000)
# ax.set_ylim(12000, 13000)
# line_center, = ax.plot([], [], color='b', label='Center', linewidth=0.5)
# # line_gc, = ax.plot([], [], color='r', label='GC', linewidth=0.75)
# ax.legend(loc='upper right')

# def animate_(i):
#     line_center.set_data(Center_df['CenterX'][i-500:i], Center_df['CenterY'][i-500:i])
#     # line_gc.set_data(GC_df['GCX'][:i], GC_df['GCY'][:i])
#     current_time = Center_df['Time'].iloc[i]
#     ax.set_title(f'Trajectory of Center at Time {current_time}')
#     return line_center,

# # Determine the number of frames for the animation
# num_frames = len(Center_df)

# ani = FuncAnimation(fig, animate_, frames=num_frames, interval=30, blit=True)

# # Save the animation
# ani.save('center_trajectory.mp4', writer='ffmpeg')
