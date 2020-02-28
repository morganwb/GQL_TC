# Run inside TC_... folder to merge procs and plot snapshots

function png2mp4(){
    cat $1 | ffmpeg \
        -y \
        -f image2pipe \
        -vcodec png \
        -r $3 \
        -i - \
        -f mp4 \
        -vcodec libx264 \
        -pix_fmt yuv420p \
        -preset slower \
        -crf 20 \
        -vf "scale=trunc(in_w/2)*2:trunc(in_h/2)*2" \
        $2
}

echo "Beginning analysis tasks..."
sleep 1
echo "Merging processor files..."
python3 -m dedalus merge_procs profiles/
python3 -m dedalus merge_procs scalar/
python3 -m dedalus merge_procs slices/
python3 -m dedalus merge_procs snapshots/
echo "Generating and animating plots..."


python3 ~/GQL_TC/python/analysis_scripts/plot_snapshots_3D.py snapshots/snapshots_s1.h5 --x_axis=theta --y_axis=z
png2mp4 'img_snapshots/*' theta_z.mp4 30
rm -r img_snapshots

python3 ~/GQL_TC/python/analysis_scripts/plot_snapshots_3D.py snapshots/snapshots_s1.h5 --x_axis=r --y_axis=z --slice=0
png2mp4 'img_snapshots/*' r_z.mp4 30
rm -r img_snapshots

python3 ~/GQL_TC/python/analysis_scripts/plot_snapshots_3D.py snapshots/snapshots_s1.h5 --x_axis=theta --y_axis=r --slice=1.5
png2mp4 'img_snapshots/*' theta_r.mp4 30
rm -r img_snapshots

echo "Analysis complete!"
