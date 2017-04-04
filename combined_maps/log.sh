    11 2017/01/30 11:22 cd -
    12 2017/01/30 11:22 ls
    13 2017/01/30 11:24 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/cube_53_56_noboundary.cm .
    14 2017/01/30 11:24 scp -r sk2534@farnam1.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/cube_53_56_noboundary.cm .
    17 2017/01/30 11:25 cd -
    18 2017/01/30 11:25 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/cube_53_56_noboundary.cm .
    10 2017/01/30 11:31 lst
    11 2017/01/30 11:31 tail log.sh
    12 2017/01/30 11:32 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/cube_53_56_noboundary.cm .
    13 2017/01/30 11:39 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/cube_53_56_noboundary.fits .
    14 2017/01/30 11:40 ls
    15 2017/01/30 11:40 which imcomb
    16 2017/01/30 11:41 imcomb in="cube_01_32_noboundary.cm,cube_33_52_noboundary.cm,cube_53_56_noboundary.cm,cube_57_90_noboundary.cm" out=cube_01_90_noboundary.cm
    17 2017/01/30 11:42 ls
    18 2017/01/30 11:42 rm *.fits
    19 2017/01/30 11:42 ls
    20 2017/01/30 11:43 fits in=cube_01_90_noboundary.cm op=xyout out=cube_01_90_noboundary.fits
    21 2017/01/30 11:43 ls
    22 2017/01/30 11:43 ds9 cube_01_90_noboundary.fits
    23 2017/01/30 11:46 lst
    24 2017/01/30 13:46 rm -rf cube_01_32_noboundary.cm/
    25 2017/01/30 13:46 rm -rf cube_33_52_noboundary.cm/
    26 2017/01/30 13:46 rm -rf cube_5*.cm
    27 2017/01/30 13:46 ls
    19 2017/01/30 16:02 sshhifi
    28 2017/02/02 9:15 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/12co.046/\*.fits .
    29 2017/02/02 9:19 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/testclean/\*.fits .
    30 2017/02/02 9:52 ls
    31 2017/02/02 9:52 ds9 12co.046.6sigma.rs.fits &
    32 2017/02/03 10:12 python
    33 2017/02/03 13:54 python
     3 2017/02/07 15:48 cd ~/Google\ Drive/12co/
     4 2017/02/07 15:48 ls
     5 2017/02/07 16:09 lst
     6 2017/02/07 16:09 du -ch cube_01_90_noboundary.cm
    10 2017/02/07 16:21 cd -
    11 2017/02/07 16:21 cp $DROPATH/python_scripts/fitscut.py .
    12 2017/02/07 16:22 ls
    13 2017/02/07 16:22 vim fitscut.py
    14 2017/02/07 16:22 vim fitscut.py
    15 2017/02/07 16:23 vim fitscut.py
    16 2017/02/07 16:23 python fitscut.py
    17 2017/02/07 16:23 vim fitscut.py
    18 2017/02/07 16:23 ds9 12co.046.6sigma.cm.fits &
    19 2017/02/07 16:24 vim fitscut.py
    20 2017/02/07 16:26 python fitscut.py
    21 2017/02/07 16:28 cp diagonal.png diagonal_6sigma.png
    22 2017/02/07 16:28 vim fitscut.py
    23 2017/02/07 16:28 python fitscut.py
    24 2017/02/07 16:28 mv diagonal.png diagonal_1sigma.png
    25 2017/02/07 16:28 ls
    26 2017/02/07 16:28 open *.png
    27 2017/02/07 16:29 vim fitscut.py
    28 2017/02/07 16:29 python fitscut.py
    29 2017/02/07 16:30 mv diagonal.png diagonal_1sigma.png
    30 2017/02/07 16:30 vim fitscut.py
    31 2017/02/07 16:30 python fitscut.py
    32 2017/02/07 16:30 mv diagonal.png diagonal_6sigma.png
    33 2017/02/07 16:31 open *.png
     3 2015/07/29 10:47 history -S
     4 2017/02/08 9:38 grep scp log.sh
     5 2017/02/08 9:39 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/testclean1sigma/12co.046.ccflux.log .
     6 2017/02/08 9:39 vim 12co.046.ccflux.log
     7 2017/02/08 9:44 cp fitscut.py plotcleanflux.py
     8 2017/02/08 9:44 vim plotcleanflux.py
     9 2017/02/08 9:44 mv plotcleanflux.py cleanflux_iteration.py
     3 2015/07/29 10:47 history -S
     4 2017/02/08 9:49 python cleanflux_iteration.py
     5 2017/02/08 9:49 python cleanflux_iteration.py
     6 2017/02/08 9:53 python cleanflux_iteration.py
     7 2017/02/08 9:54 python cleanflux_iteration.py
     8 2017/02/08 9:59 ls *.png
     9 2017/02/08 10:00 python cleanflux_iteration.py
    10 2017/02/08 9:45 vim cleanflux_iteration.py
    11 2017/02/08 10:02 mv cleanflux_iteration.py plotcleanflux.py
    10 2017/02/08 10:02 ls
    11 2017/02/08 10:02 ls *.py
    12 2017/02/08 10:02 mv *.py $DROPATH/python_scripts/
    34 2017/02/08 10:14 ls
    35 2017/02/08 10:14 ds9 12co.046.6sigma.cm.fits
    35 2017/02/08 10:14 ds9 12co.046.6sigma.cm.fits
    36 2017/02/08 10:18 open .
    37 2017/02/08 10:21 cp $DROPATH/python_scripts/plotcleanflux.py .
    38 2017/02/08 10:21 vim plotcleanflux.py
    39 2017/02/08 10:21 python plotcleanflux.py
    40 2017/02/08 10:24 mv plotcleanflux.py $DROPATH/python_scripts/
    13 2017/02/08 11:02 ls
    14 2017/02/08 11:53 ds9 12co.046.1sigma.rs.fits
    15 2017/02/08 11:53 which smir
    16 2017/02/08 11:53 source ~/miriad_tcsh
    17 2017/02/08 11:53 ds9 12co.046.1sigma.rs.fits
    66 2017/02/09 11:49 cd ~/Google\ Drive/12co/
    67 2017/02/09 11:49 ls
    68 2017/02/09 11:49 ds9 new_cube_01_90_noboundary.fits &
    78 2017/02/11 20:57 cd ~/Google\ Drive/12co/
    79 2017/02/11 20:57 ls
    80 2017/02/11 20:57 ds9 new_cube_01_90_noboundary.fits &
    81 2017/02/11 23:03 python
    82 2017/02/12 11:49 grep scp log.sh
    83 2017/02/12 11:50 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/test058/12co.058.cm.fits .
    12 2017/02/12 11:50 tail log.sh
    13 2017/02/12 11:50 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/test058/12co.058.cm.fits .
    14 2017/02/12 11:52 which ds9
    15 2017/02/12 11:52 smir
    16 2017/02/12 11:52 source ~/miriad_tcsh
    17 2017/02/12 11:52 ds9 12co.058.cm.fits
    18 2017/02/12 11:58 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co/test058/12co.058.map.fits .
    19 2017/02/12 12:00 ds9 12co.058.map.fits
    20 2017/02/12 12:02 ls
    84 2017/02/12 12:03 cp $research/2016/NRO/12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits .
    85 2017/02/12 12:03 mv $research/2016/NRO/12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits .
    86 2017/02/12 12:03 rm -rf $research/2016/NRO/12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.mir/
    87 2017/02/12 12:03 ls
    88 2017/02/12 12:04 ds9 12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits
    89 2017/02/12 12:05 mv $research/2016/NRO/12CO_*.fits .
    90 2017/02/12 12:05 ds9 12CO_specsmooth_nonan.fits
    91 2017/02/12 12:06 ds9 new_cube_01_90_noboundary.fits &
    91 2017/02/12 12:06 ds9 new_cube_01_90_noboundary.fits &
    92 2017/02/12 12:07 ds9 12CO_specsmooth_nonan.fits
    93 2017/02/12 14:35 ds9 new_cube_01_90_noboundary.fits &
    94 2017/02/12 17:07 exit
     1 2017/02/12 18:45 cd Google\ Drive/12co/
     2 2017/02/12 18:45 ls
     3 2017/02/12 18:45 grep scp log.sh
     4 2017/02/12 18:45 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
     5 2017/02/12 19:02 ds9 12co.002.map.fits &
     6 2017/02/12 19:15 ds9 12co.002.map.fits &
     7 2017/02/12 22:50 ds9 12co.002.map.fits &
     8 2017/02/12 23:47 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
     8 2017/02/12 23:47 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
     9 2017/02/12 23:48 ds9 12co.002.map.fits &
     9 2017/02/12 23:48 ds9 12co.002.map.fits &
    10 2017/02/13 8:02 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    10 2017/02/13 8:02 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    11 2017/02/13 8:02 ds9 12co.002.map.fits &
    11 2017/02/13 8:02 ds9 12co.002.map.fits &
    12 2017/02/13 9:24 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    13 2017/02/13 9:25 ds9 12co.002.map.fits &
    14 2017/02/13 10:36 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    15 2017/02/13 10:37 ds9 12co.002.map.fits &
    16 2017/02/13 11:53 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    17 2017/02/13 11:54 ds9 12co.002.map.fits &
    22 2017/02/13 11:56 sshfarnam
    18 2017/02/13 12:58 grep grace log.sh
    19 2017/02/13 13:00 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testsmallimsize/12co.002/12co.002.cm.fits .
    23 2017/02/13 13:01 tail log.sh
    24 2017/02/13 13:01 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testsmallimsize/12co.002/12co.002.cm.fits .
    20 2017/02/13 13:22 lst
    21 2017/02/13 13:23 ds9 12co.002.cm.fits
    22 2017/02/14 8:26 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testsmallimsize/12co.002/12co.002.cm.fits .
    25 2017/02/14 8:26 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testsmallimsize/12co.002/12co.002.cm.fits .
    23 2017/02/14 8:58 lst
    24 2017/02/14 8:58 ds9 12co.002.cm.fits
    25 2017/02/14 16:31 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    26 2017/02/14 16:32 ds9 12co.002.cm.fits
    27 2017/02/14 16:33 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    28 2017/02/14 16:34 ds9 12co.002.cm.fits
    29 2017/02/14 16:35 ls
    30 2017/02/14 16:35 lst
    31 2017/02/14 16:35 rm 12co.002.*
    32 2017/02/14 16:35 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    33 2017/02/14 16:36 lst
    34 2017/02/14 16:37 ds9 12co.002.map.fits
    35 2017/02/14 21:01 scp carmaorion@hifi.caltech.edu:/hifi/carmaorion/orion/images/sk/12co/12co.002/12co.002.map.fits .
    36 2017/02/14 21:04 ds9 12co.002.map.fits
    26 2017/02/14 21:10 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize100/12co.002/12co.002.cm.fits .
    27 2017/02/14 21:13 ds9 12co.002.cm.fits
    28 2017/02/14 21:17 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize100/12co.002/12co.002.map.fits .
    29 2017/02/14 21:20 ds9 12co.002.map.fits
    30 2017/02/14 21:21 ds9 12co.002.cm.fits
    30 2017/02/14 21:21 ds9 12co.002.cm.fits
    31 2017/02/15 11:34 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120/12co.002/12co.002.map.fits .
    32 2017/02/15 11:34 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120/12co.002/12co.002.cm.fits .
    33 2017/02/15 11:38 ds9 12co.002.cm.fits
    34 2017/02/15 21:37 ds9 12co.002.cm.fits
    35 2017/02/15 21:39 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120robust-2/12co.002/12co.002.cm.fits .
    36 2017/02/15 21:41 ds9 12co.002.cm.fits
    37 2017/02/16 7:30 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120robust0.5sample3500/12co.002/12co.002.cm.fits .
    38 2017/02/16 7:31 ds9 12co.002.cm.fits
    39 2017/02/16 7:31 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120robust0.5sample2300/12co.002/12co.002.cm.fits .
    40 2017/02/16 7:33 ds9 12co.002.cm.fits
    41 2017/02/16 7:34 scp sk2534@grace.hpc.yale.edu:/project/fas/arce/sk2534/12co/testimsize120robust-2sample3500/12co.002/12co.002.cm.fits .
    42 2017/02/16 7:35 ds9 12co.002.cm.fits
    43 2017/02/16 16:41 ls
    44 2017/02/16 17:31 exit
    34 2017/02/05 10:24 python
     4 2017/02/20 9:15 cd ..
     5 2017/02/20 9:15 ds9 12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits
     6 2017/02/20 9:15 source ~/miriad_tcsh
     7 2017/02/20 9:15 ds9 12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits
     8 2017/02/20 9:16 python
     9 2017/02/24 13:59 ls
    10 2017/02/24 13:59 ds9 new_cube_01_90_noboundary.fits
    11 2017/02/24 14:06 ds9 new_cube_01_90_noboundary.fits
    12 2017/02/24 16:48 ds9 new_cube_01_90_noboundary.fits
    13 2017/02/25 10:13 mkdir test046imsize
    14 2017/02/25 10:13 ls
     4 2017/02/25 10:24 cd ../../12co/
     5 2017/02/25 10:24 ls
     4 2017/02/28 9:00 cd ~/Google\ Drive/12co/
     5 2017/02/28 9:00 grep scp log.sh
     6 2017/02/28 9:00 grep scp test046imsize/log.sh
     7 2017/02/28 9:00 scp -r sk2534@farnam.hpc.yale.edu:/ysm-gpfs/project/sk2534/12co_31_60/12co.006/12co.036.cm.fits .
     8 2017/02/28 9:03 which ds9
     9 2017/02/28 9:04 source ~/miriad_tcsh
    10 2017/02/28 9:04 ds9 12co.036.cm.fits
     4 2017/03/01 9:24 cd ~/Google\ Drive/12co/
     5 2017/03/01 9:24 ls
     6 2017/03/01 9:24 ds9 12co_31_60.cm.fits
     7 2017/03/01 9:40 sshkong
     8 2017/03/01 15:24 grep gauss $DROPATH/python_scripts/*.py
     9 2017/03/01 15:25 python
    10 2017/03/01 18:41 sshhifi
     3 2015/07/29 10:47 history -S
    11 2017/03/02 11:18 sshkong
    12 2017/03/03 9:53 ds9 12co_pix_2.cm.fits
    13 2017/03/03 9:54 source ~/miriad_tcsh
    14 2017/03/03 9:54 ds9 12co_pix_2.cm.fits
    15 2017/03/03 9:55 ls
    12 2017/03/03 9:56 cd ~/GoogleDrive/12co/
    13 2017/03/03 9:56 ls
    16 2017/03/03 9:58 python
    17 2017/03/03 10:24 python
    14 2017/03/03 9:57 ds9 12co_pix_1p5.cm.fits
     3 2015/07/29 10:47 history -S
    15 2017/03/04 11:39 python
    47 2017/03/14 12:34 cd ../../12co/
     3 2015/07/29 10:47 history -S
     4 2017/03/14 12:37 open .
    48 2017/03/14 12:34 ds9 12co_pix_2.cm.fits
     5 2017/03/18 14:23 grep scp log.sh
     6 2017/03/18 14:24 scp sk2534@grace-next.hpc.yale.edu:/project/fas/arce/sk2534/13co/13co.001/13co.001.cm.fits .
     7 2017/03/18 14:25 scp sk2534@grace-next.hpc.yale.edu:/project/fas/arce/sk2534/13co/13co/13co.001/13co.001.cm.fits .
     8 2017/03/18 14:36 source ~/miriad_tcsh
     9 2017/03/18 14:36 ds9 13co.001.cm.fits
    10 2017/03/18 20:11 scp sk2534@grace-next.hpc.yale.edu:/project/fas/arce/sk2534/13co/13co/13co.001/13co.001.cm.fits .
    11 2017/03/18 20:16 ds9 13co.001.cm.fits
    12 2017/03/21 14:48 ds9 13co.001.cm.fits
    12 2017/03/21 14:48 ds9 13co.001.cm.fits
     4 2017/03/22 16:25 cd ~/Google\ Drive/12co/
     5 2017/03/22 16:25 ls
     6 2017/03/22 16:25 lst
     7 2017/03/22 16:26 ds9 13co.001.cm.fits
     8 2017/03/22 16:26 source ~/miriad_tcsh
     9 2017/03/22 16:26 ds9 13co.001.cm.fits
     2 2017/03/24 9:12 cd GoogleDrive/12co/
     3 2017/03/24 9:12 ls
     4 2017/03/24 9:12 ds9 12CO_20161017_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms_YS.fits
     5 2017/03/24 9:27 sshhifi
     6 2017/03/24 9:27 ds9 12CO_specsmooth_nonan.fits
    11 2017/03/24 10:39 cd ../12co/
    12 2017/03/24 10:39 ls
    13 2017/03/24 10:39 ds9 12co_pix_2.cm.fits
    14 2017/03/24 10:58 ds9 12co_pix_2.cm.fits
     4 2017/03/27 20:42 cd ~/Google\ Drive/12co/
     5 2017/03/27 20:42 l
     6 2017/03/27 20:42 ls
     7 2017/03/27 20:42 rm -rf test046imsize/
     8 2017/03/27 20:42 rm -rf test058/
     9 2017/03/27 20:43 rm 12co.0*
    10 2017/03/27 20:43 ls
    11 2017/03/27 20:43 rm 13co.001.cm.fits
    12 2017/03/27 20:51 ls
    13 2017/03/27 20:51 rm Icon
    14 2017/03/27 20:51 ls
    15 2017/03/27 22:02 ds9 12co_pix_2.cm.fits
    16 2017/03/27 22:02 source ~/miriad_tcsh
    17 2017/03/27 22:02 ds9 12co_pix_2.cm.fits
    17 2017/03/27 22:02 ds9 12co_pix_2.cm.fits
    20 2017/03/27 23:39 cd -
    21 2017/03/27 23:39 ds9 12co_pix_2.cm.fits
    22 2017/03/28 0:20 which convert
    23 2017/03/28 0:20 convert
    24 2017/03/28 0:20 man convert
    25 2017/03/28 0:21 ls
     3 2015/07/29 10:47 history -S
    26 2017/03/28 0:57 which convert
     7 2017/03/28 0:57 cd ~/Google\ Drive/12co/
     8 2017/03/28 0:57 ls
     2 2017/03/28 1:06 cd ~/GoogleDrive/12co/
     3 2017/03/28 1:06 ls
     4 2017/03/28 1:06 rm -rf carmaorion12co.mp4
    30 2017/03/28 8:34 cd ..
    31 2017/03/28 8:34 ds9 12co_pix_2.cm.fits
    32 2017/03/28 8:34 source ~/miriad_tcsh
    33 2017/03/28 8:34 ds9 12co_pix_2.cm.fits
    33 2017/03/28 8:34 ds9 12co_pix_2.cm.fits
    34 2017/03/28 10:50 cd ../13
     7 2017/03/28 12:03 cd ../
     8 2017/03/28 12:03 ls
    11 2017/03/28 12:04 cd ..
    12 2017/03/28 12:04 mkdir channels_13co_pix_2/
    13 2017/03/28 16:31 mv channels_13co_pix_2 ../13co/
