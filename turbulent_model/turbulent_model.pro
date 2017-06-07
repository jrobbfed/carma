pro turbulent_model, outfile=outfile, dist=dist,$
                     pix_size=pix_size, vstep=vstep,$
                     acen=acen, dcen=dcen, thickness=thickness,$
                     fwhm=fwhm, beta=beta, r=r, dr=dr,$
                     vexp=vexp, depth_offset=depth_offset,$
                     vel_offset=vel_offset, v0=v0,$
                     ignore_cloud=ignore_cloud,$
                     write_fits=write_fits,$
                     save_ppv=save_ppv

  ;- EDIT THESE NUMBERS
;   dist = 414.                       ;- distance to Orion A, pc -  Menten + 2007
;   pix_size = 7.5                ;- nyquist pixel scale, arcsec - NRO maps
;   vstep = 0.099                 ;- velocity resolution, km/s - 12CO NRO 
;   acen = 83.72707               ;- RA at bubble center, shell 18
;   dcen = -5.07792               ;- Dec at bubble center, shell 18

;   thickness = 1.0                 ;- cloud thickness, pc, from ~size of largest shell.
;   fwhm = 1.0                  ;- cloud velocity fwhm, km/s
;   beta = 2.0                    ;- spectral index of velocity field
  
; ;  r = 0.31                       ;- radius of bubble, pc
;   r = 0.31                       ;- radius of bubble, pc
;   dr = 0.2                     ;- thickness of bubble, pc
;   vexp = 1.6                     ;- km/s -- cap-to-midplane
; ;  vexp = 3.                     ;- km/s -- cap-to-midplane
;   depth_offset = 0.0            ;- offset of bubble from cloud center, in pc
;   vel_offset = 0.0              ;- difference between cloud and bubble velocity, km/s
;   v0 = 13.6;- systematic velocity (mid-channel of output cube)

;  if N_elements(outdir) then outdir = './'
  if N_elements(outfile) eq 0 then outfile = 'turb_ppv.fits'
  if N_elements(dist) eq 0 then dist = 414.
  if N_elements(pix_size) eq 0 then pix_size = 7.5
  if N_elements(vstep) eq 0 then vstep = 0.099
  if N_elements(acen) eq 0 then acen = 83.72707
  if N_elements(dcen) eq 0 then dcen = -5.07792
  if N_elements(thickness) eq 0 then thickness = 1.0
  if N_elements(fwhm) eq 0 then fwhm = 4.0 
  if N_elements(beta) eq 0 then beta = 0.0
  if N_elements(r) eq 0 then r = 0.31
  if N_elements(dr) eq 0 then dr = 0.2
  if N_elements(vexp) eq 0 then vexp = 1.6
  if N_elements(depth_offset) eq 0 then depth_offset = 0.0 
  if N_elements(vel_offset) eq 0 then vel_offset = 0.0
  if N_elements(v0) eq 0 then v0 = 13.6
  if N_elements(ignore_cloud) eq 0 then ignore_cloud = 0
  if N_elements(write_fits) eq 0 then write_fits = 1
  if N_elements(save_ppv) eq 0 then save_ppv = 0


  ;- do computation at pixel scale 2x finer than end result
  scale = pix_size / 206265. * dist / 2 ;- pc per pixel
  print, scale, " pc/pixel"

  ;- dimensions of working ppp cubes
  if ignore_cloud then begin
    sz = floor([4 * r / scale, 4. * r / scale, 4. * r / scale])
  endif else begin
    sz = floor([4 * r / scale, 4. * r / scale, 1.1 * thickness / scale])
  endelse

  den = fltarr(sz[0], sz[1], sz[2])
  indices, den, x, y, z, center = fltarr(3)

  rr = sqrt(x^2 + y^2 + (z - depth_offset / scale)^2) * scale

  if ignore_cloud then begin
    in_slab = intarr(sz[0], sz[1], sz[2]) + 1 ; All ones.
  endif else begin 
    in_slab = abs(z) * scale lt (thickness / 2.)
  endelse

  in = rr lt (r - dr/2.) and in_slab
  on = rr ge (r - dr/2.) and rr lt (r + dr/2.) and in_slab
  out = rr ge (r + dr/2.) and in_slab

  print, total(in), total(rr lt r - dr/2.)
  ;print, minmax(rr)
  ;print, minmax(z), depth_offset / scale

  ;- density field: shell expanding into uniform exterior
  ;- evacuated material uniformly redistributed over shell boundary
  ;- units are arbitrary
  print, "Computing Density of Shell"
  

  ratio = total(in) / total(on)

  if ignore_cloud then begin
    den = on
  endif else begin
    den = out + on * (1 + ratio)
  endelse
  print, "Generating turbulent cloud velocity field"

  ;velocity field of bubble
  vel = on * (vexp * (z * scale - depth_offset) / rr + vel_offset)
  bad = where(rr eq 0, bad_ct)
  if bad_ct ne 0 then vel[bad] = 0.

  if ~ignore_cloud then begin
    ;- superpose turbulent velocity field of cloud. Units are km/s
    seed = 101
    vel = vel + cloud(sz, beta = beta, seed = seed) * fwhm / 2.35
  endif

  ;- velocity spacing of output cube
  vlo = -3 * fwhm < vel_offset - 1.2 * vexp
  vhi = 3 * fwhm > vel_offset + 1.2 * vexp
  vcen = arrgen(vlo, vhi, vstep)

  ;- ppv cube
  print, "Gridding PPV cube."

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
  
  if save_ppv then begin
    save, ppv, outfile
  endif

  if write_fits then begin
  ;- make a header
    print, "Writing fits."
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
    sxaddpar, hdr, 'CRVAL3', v0, 'KM/S'
    sxaddpar, hdr, 'CDELT3', vstep, 'KM/S'
    sxaddpar, hdr, 'CUNIT3', 'km/s'

    sxaddpar, hdr, 'THICK', thickness, 'Cloud Thickness (pc)'
    sxaddpar, hdr, 'DIST', dist, 'Distance to cloud (pc)'
    sxaddpar, hdr, 'V_FWHM', fwhm, 'Cloud velocity spread (km/s)'
    sxaddpar, hdr, 'BETA', beta, 'Cloud velocity power spectrum index'
    sxaddpar, hdr, 'VEXP', vexp, 'Expansion vel (km/s - cap to midplane)'
    sxaddpar, hdr, 'R', r, 'Bubble size (pc)'
    sxaddpar, hdr, 'DR', dr, 'Bubble thickness (pc)'
    sxaddpar, hdr, 'ZOFF', depth_offset, 'Bubble depth offset (pc)'
    sxaddpar, hdr, 'VOFF', vel_offset, 'Bubble-Cloud vel (km/s)'

    ;- write out result
    writefits, outfile, ppv, hdr ;- the ppv cube
  endif
  writefits, "den.fits", den
  writefits, "vel.fits", vel

;- write out ppv, den, vel as a multi-extension cube
;  file = 'all_'+outfile
;  mwrfits, ppv, file, hdr, /create
;
;  ;- write out density, velocity cubes
;  mkhdr, hdr, den
;  sz = size(den)
;  sxaddpar, hdr, 'CTYPE1', 'X----CAR'
;  sxaddpar, hdr, 'CRPIX1', sz[1]/2.
;  sxaddpar, hdr, 'CRVAL1', 0, 'PC'
;  sxaddpar, hdr, 'CDELT1', scale, 'PC'
;
;  sxaddpar, hdr, 'CTYPE2', 'Y----CAR'
;  sxaddpar, hdr, 'CRPIX2', sz[2]/2.
;  sxaddpar, hdr, 'CRVAL2', 0, 'PC'
;  sxaddpar, hdr, 'CDELT2', scale, 'PC'
;
;  sxaddpar, hdr, 'CTYPE3', 'Z----CAR'
;  sxaddpar, hdr, 'CRPIX3', sz[3]/2.
;  sxaddpar, hdr, 'CRVAL3', 0, 'PC'
;  sxaddpar, hdr, 'CDELT3', scale, 'PC'
;
;  sxaddpar, hdr, 'BUNIT', 'cm^-3', 'DENSITY'
;
;  mwrfits, den, file, hdr ;- a higher-res 3D density cube
;
;  sxaddpar, hdr, 'BUNIT', 'km/s', 'RADIAL VELOCITY'
;  mwrfits, vel, file, hdr ;- higher-res 3D radial velocity cube
end
