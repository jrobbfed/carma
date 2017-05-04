pro turbulent_model, outfile = outfile

  ;- EDIT THESE NUMBERS
  d = 300                       ;- distance to cloud, pc
  pix_size = 23                 ;- nyquist pixel scale, arcsec
  vstep = 0.1                   ;- velocity resolution, km/s
  acen = 52.125                 ;- RA at bubble center
  dcen = 31.3                   ;- Dec at bubble center

  thickness = 1.2               ;- cloud thickness, pc
  fwhm = 1.6                    ;- cloud velocity fwhm, km/s
  beta = 2.0                    ;- spectral index of velocity field
  
  r = 0.6                       ;- radius of bubble, pc
  dr = 0.25                     ;- thickness of bubble, pc
  vexp = 3.                     ;- km/s -- cap-to-midplane
  depth_offset = 0.0            ;- offset of bubble from cloud center, in pc
  vel_offset = 0.0              ;- difference between cloud and bubble velocity, km/s

  ;- STOP EDITING
  if ~keyword_set(outfile) then outfile = 'turb_ppv.fits'


  ;- do computation at pixel scale 2x finer than end result
  scale = pix_size / 206265. * d / 2 ;- pc per pixel
  ;- dimensions of working ppp cubes
  sz = floor([4 * r / scale, 4. * r / scale, 1.1 * thickness / scale])

  den = fltarr(sz[0], sz[1], sz[2])
  indices, den, x, y, z, center = fltarr(3)

  rr = sqrt(x^2 + y^2 + (z - depth_offset / scale)^2) * scale
  in_slab = abs(z) * scale lt (thickness / 2)
  in = rr lt (r - dr/2.) and in_slab
  on = rr ge (r - dr/2.) and rr lt (r + dr/2.) and in_slab
  out = rr ge (r + dr/2.) and in_slab

  print, total(in), total(rr lt r - dr/2.)
  print, minmax(rr)
  print, minmax(z), depth_offset / scale

  ;- density field: shell expanding into uniform exterior
  ;- evacuated material uniformly redistributed over shell boundary
  ;- units are arbitrary
  ratio = total(in) / total(on)
  den = out + on * (1 + ratio)

  ;- turbulent velocity field of cloud. Units are km/s
  seed = 100
  vel = cloud(sz, beta = beta, seed = seed) * fwhm / 2.35
  ;- superpose velocity field of bubble
  vel_bubble = on * (vexp * (z * scale - depth_offset) / rr + vel_offset)
  bad = where(rr eq 0, bad_ct)
  if bad_ct ne 0 then vel_bubble[bad] = 0.
  vel += vel_bubble


  ;- velocity spacing of output cube
  vlo = -3 * fwhm < vel_offset - 1.2 * vexp
  vhi = 3 * fwhm > vel_offset + 1.2 * vexp
  vcen = arrgen(vlo, vhi, vstep)

  ;- ppv cube
  ppv = ppp2ppv(den, vel, vcen)

  ;- downsample by 2x to requested pixel scale
  sz = size(ppv)
  ppv = congrid(ppv, sz[1]/2, sz[2]/2, sz[3])

  ;- convolve with a psf of fwhm = 2 pixels
  psf = psf_gaussian(npixel = 10, fwhm = 2., /norm)
  psf /= total(psf)
  sz = size(ppv)
  for i = 0, sz[3]-1, 1 do begin
     ppv[*,*,i] = convolve(reform(ppv[*,*,i]), psf, ft_psf = ft_psf)
  endfor

  ;- make a header
  mkhdr, hdr, ppv

  sxaddpar, hdr, 'CTYPE1', 'RA---TAN'
  sxaddpar, hdr, 'CRPIX1', sz[1]/2.
  sxaddpar, hdr, 'CRVAL1', acen, 'DEGREES'
  sxaddpar, hdr, 'CDELT1', pix_size / 3600., 'DEGREES'

  sxaddpar, hdr, 'CTYPE2', 'DEC--TAN'
  sxaddpar, hdr, 'CRPIX2', sz[2]/2.
  sxaddpar, hdr, 'CRVAL2', dcen, 'DEGREES'
  sxaddpar, hdr, 'CDELT2', pix_size / 3600., 'DEGREES'

  sxaddpar, hdr, 'CTYPE3', 'VELO-LSR'
  sxaddpar, hdr, 'CRPIX3', sz[3]/2.
  sxaddpar, hdr, 'CRVAL3', 0.0, 'KM/S'
  sxaddpar, hdr, 'CDELT3', vstep, 'KM/S'

  sxaddpar, hdr, 'THICK', thickness, 'Cloud Thickness (pc)'
  sxaddpar, hdr, 'DIST', d, 'Distance to cloud (pc)'
  sxaddpar, hdr, 'V_FWHM', fwhm, 'Cloud velocity spread (km/s)'
  sxaddpar, hdr, 'BETA', beta, 'Cloud velocity power spectrum index'
  sxaddpar, hdr, 'VEXP', vexp, 'Expansion vel (km/s - cap to midplane)'
  sxaddpar, hdr, 'R', r, 'Bubble size (pc)'
  sxaddpar, hdr, 'DR', dr, 'Bubble thickness (pc)'
  sxaddpar, hdr, 'ZOFF', depth_offset, 'Bubble depth offset (pc)'
  sxaddpar, hdr, 'VOFF', vel_offset, 'Bubble-Cloud vel (km/s)'

  ;- write out result
  writefits, outfile, ppv, hdr ;- the ppv cube


  ;- write out ppv, den, vel as a multi-extension cube
  file = 'all_'+outfile
  mwrfits, ppv, file, hdr, /create

  ;- write out density, velocity cubes
  mkhdr, hdr, den
  sz = size(den)
  sxaddpar, hdr, 'CTYPE1', 'X----CAR'
  sxaddpar, hdr, 'CRPIX1', sz[1]/2.
  sxaddpar, hdr, 'CRVAL1', 0, 'PC'
  sxaddpar, hdr, 'CDELT1', scale, 'PC'

  sxaddpar, hdr, 'CTYPE2', 'Y----CAR'
  sxaddpar, hdr, 'CRPIX2', sz[2]/2.
  sxaddpar, hdr, 'CRVAL2', 0, 'PC'
  sxaddpar, hdr, 'CDELT2', scale, 'PC'

  sxaddpar, hdr, 'CTYPE3', 'Z----CAR'
  sxaddpar, hdr, 'CRPIX3', sz[3]/2.
  sxaddpar, hdr, 'CRVAL3', 0, 'PC'
  sxaddpar, hdr, 'CDELT3', scale, 'PC'

  sxaddpar, hdr, 'BUNIT', 'cm^-3', 'DENSITY'

  mwrfits, den, file, hdr ;- a higher-res 3D density cube

  sxaddpar, hdr, 'BUNIT', 'km/s', 'RADIAL VELOCITY'
  mwrfits, vel, file, hdr ;- higher-res 3D radial velocity cube
end
