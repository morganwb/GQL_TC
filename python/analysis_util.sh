# Run inside TC_... folder to merge procs and plot snapshots

echo "Beginning analysis tasks..."
sleep 1
echo "Merging processor files..."
python3 -m dedalus merge_procs profiles/
python3 -m dedalus merge_procs scalar/
python3 -m dedalus merge_procs slices/
python3 -m dedalus merge_procs snapshots/
echo "Generating plots..."
python3 ~/GQL_TC/python/plot_snapshots_3D.py snapshots/snapshots_s1.h5
#echo "Animating plots..."
#ffmpeg "img_snapshots/*" frames.mp4 30
echo "Analysis complete!"
