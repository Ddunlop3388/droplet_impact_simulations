qcc -Wall -O2 -grid=octree droplet_impact.c -o di -L$BASILISK/gl -lglutils -lfb_tiny -lm

./di
