function tile {
    gmic $1.png $1.png $1.jpg $1_restored.png -remove_opacity[0-3] -crop[0] 0,0,512,512 -crop[1-3] $2,$3,$(expr $2 + 64),$(expr $3 + 64) _interpolation=0 -resize[1-3] 512,512 -append_tiles[0-3] , -o $1_tiles_.png
    pngcrush $1_tiles_.png $1_tiles.png
    rm $1_tiles_.png
}

tile lena 230 240
tile deviantart 230 150
