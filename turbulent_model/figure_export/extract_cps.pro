;+ 
; Extract a postage stamp image of a CPS from 
; the big fits images
;-
pro extract_cps

  co_12_file = '../../PerA_12coFCRAO_F_xyv_TAN.fits'
  co_13_file = '../../PerA_13coFCRAO_F_xyv_TAN.fits'

  ;- offset of image center from CPS center, units of radii.
  ;- one for each CPS
  off_x = [0., 0., 0., 0., $
           0., 0., 0., 0., $
           0., 0., 0., 0.]
  off_y = [0., 0., 0., 0., $
           0., 0., 0., 0., $
           0., 0., 0., 0]

  ;- size of postage stamp, in CPS radii
  wid = [4, 4., 4., 4., $
         4., 4., 4., 4., $
         4., 4., 4., 4.]


  ;- positions of cps
  cps = cps()
  ra = cps.ra
  dec = cps.dec
  radius = cps.r


  if ~file_test(co_12_file) then $
     message, 'Could not find 12CO data cube: ' + co_12_file
  if ~file_test(co_13_file) then $
     message, 'Could not find 13CO data cube: ' + co_13_file

  common extract_common, m, h, m12, h12
  if n_elements(h12) eq 0 then begin
     m = mrdfits(co_13_file,0,h)
     m12 = mrdfits(co_12_file,0,h12)
  endif

  sz = size(m)
                                ;for i = 0, n_elements(ra) - 1, 1 do begin
  for i = 0, 11, 1 do begin
     print, i

     acen = ra[i] + off_x[i] * radius[i] / cos(dec[i] * !dtor)
     dcen = dec[i] + off_y[i] * radius[i]

     cdelt = abs(sxpar(h, 'cdelt1'))
     fov = wid[i] * radius[i]
     npix = ceil(fov / cdelt)
     
     zlo = 103
     zhi = 350
     zstep = zhi - zlo + 1
     image = fltarr(npix, npix, zstep)
     mkhdr, h2, image
     sxaddpar, h2, 'CTYPE1', 'RA---CAR'
     sxaddpar, h2, 'CTYPE2', 'DEC--CAR'
     sxaddpar, h2, 'CTYPE3', 'VELO-LSR'
     sxaddpar, h2, 'crval1', acen
     sxaddpar, h2, 'crpix1', (npix+1)/2.
     sxaddpar, h2, 'cdelt1', -cdelt
     sxaddpar, h2, 'crval2', dcen
     sxaddpar, h2, 'crpix2', (npix+1) / 2.
     sxaddpar, h2, 'cdelt2', cdelt
     sxaddpar, h2, 'crval3', sxpar(h, 'crval3')
     sxaddpar, h2, 'crpix3', sxpar(h, 'crpix3') - zlo
     sxaddpar, h2, 'cdelt3', sxpar(h, 'cdelt3')


     indices, image, xx, yy, zz
     xyzadv, h2, xx, yy, zz, aa, dd, vv, /ortho
     advxyz, h, aa, dd, vv, xx, yy, zz, /ortho
     image = interpolate(m, xx, yy, zz)
     
     writefits, 'cps_13co_'+string(i, format='(i2.2)')+'.fits', image, h2

     advxyz, h12, aa, dd, vv, xx, yy, zz, /ortho
     image = interpolate(m12, xx, yy, zz)
     writefits, 'cps_12co_'+string(i, format='(i2.2)')+'.fits', image, h2
  endfor
     
end
