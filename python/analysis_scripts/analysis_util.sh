# Run inside TC_... folder to merge procs and plot snapshots

echo "Beginning analysis tasks..."
sleep 1
echo "Merging processor files..."
python3 -m dedalus merge_procs profiles/
python3 -m dedalus merge_procs scalar/
python3 -m dedalus merge_procs slices/
python3 -m dedalus merge_procs snapshots/
echo "Generating plots..."
python3 ~/GQL_TC/python/analysis_scripts/plot_snapshots_3D.py snapshots/snapshots_s1.h5
echo "Animating plots..."
ffmpeg -i img_snapshots/* -y -f image2pipe -vcodec png -r 30 -f mp4 -vcodec libx264 -pix_fmt yuv420p -preset slower -crf 10 -vf "scale=trunc(in_w/2)*2:trunc(in_h/2)*2" "snapshots.mp4"
echo "Analysis complete!"
