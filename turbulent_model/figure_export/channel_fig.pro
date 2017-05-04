;+
; PURPOSE:
;  Standardized routine to make channel maps for a CPS source. The
;  plots have the following characteristics:
;
;  - Each shows 9 velocities, evenly spaced from the front to the back
;    of the shell
;  - The candidate driving source is plotted as a blue star
;  - The profile of an expanding sphere is overplotted as a red
;    annulus
;  - The field of view is a constant number of bubble radii (about 4)
;
; If you want to use 13CO data instead of 12CO, call channel_fig,
; /co13
;-
pro channel_fig, number, co13 = co13

  ;- contour levels for each bubble
  contour_start = [0.4, 0.4, 6, 0.75, $
                   0.6, 1.5, 1.0, 1.0, $
                   0.8, 0.8, 0.5, 0.9]
  contour_step = [0.4, 0.4, 1, 0.75 * 2, $
                  1.2 * 2, 1.5, 1.0 * 2, 2, $
                  1.6 ,1.6, 1.0, 0.9]
  lev = contour_start[number] + findgen(20) * contour_step[number]

  ;- location of driving source candidates
  src_ra = [ten(3, 27, 40), -1, ten(3, 29, 19.8), -1, $
            ten(3, 39, 44.7), ten(3, 41, 22.1), ten(3, 42, 10.9), ten(3, 44, 19.1), $
            ten(3 ,44, 19.2), ten(3, 44, 34.2), ten(3, 44, 50.7), ten(3, 47, 47.1)]
  src_dec = [ten(31, 15, 40), -1, ten(31, 24, 57), -1, $
             ten(31, 55, 33), ten(31, 54, 34), ten(32, 08, 17), ten(32, 17, 18), $
             ten(32, 07, 35), ten(32, 09, 46), ten(32, 19, 06), ten(33, 04, 04)]
  
  ;- load data
  fmt = keyword_set(co13) ? '("cps_13co_", i2.2, ".fits")' : '("cps_12co_", i2.2, ".fits")'
  co = string(number, format= fmt)

  co = mrdfits(co, 0, h)

  ;- a circle
  theta = arrgen(0., 2 * !pi, nstep = 50)
  xcircle = cos(theta)
  ycircle = sin(theta)

  ;- get ra, dec, vcen, vexp, r, dr for the bubble
  p = cps(number)
  pixelize, p, h ;- convert to pixel coordinates

  ;- set greyscale
  nanswap, co, 0
  co_raw = co
  floor = where(co lt contour_start[number])
  co = sigrange(co, range = .98) > contour_start[number]
  co -= contour_start[number]
  co /= max(co)
  co = 1 - co
  co = byte(co * 150) + (253B - 150B)
  co[floor] = 254B

  ;- set up plotting grid
  nrow = 3
  ncol = 3
  ntile = nrow * ncol
  xborder = .02
  yborder = .02
  xoff = .10
  yoff = .10
  sz = size(co)
  aspect_ratio = 1. * sz[2] / sz[1]
  dx = (1. - xoff) / ncol
  dy = (1. - yoff) / ncol * aspect_ratio

  ;- set up postscript device
  set_plot, 'ps'
  loadct, 0, /silent
  device, file='channel_'+string(number+1, format='(i2.2)')+'.eps', $
          /encap, /color, bits_per_pixel = 8, $
          xsize = 7, ysize = 7, /in, $
          /helvetica
  !p.font = 0
  !p.charsize = .9
  tvlct, 220, 20, 60, 255
  tvlct, 65, 105, 225, 254

  ;- define the "star" symbol to use when plotting driving candidate
  xstar = [0, .225, .951, .365, .5877, 0, -.5877, -.36, -.951056, -0.22, 0]
  ystar = [1, .30917, .30917, -.120, -0.809, -.384, -.809017, -.120, .30917, .30917, 1]
  usersym, xstar, ystar, /fill
     
  ;- generate each velocity slice
  ;- 9 velocity slices, stepping evenly from front to back of shell
  nslice = nrow * ncol
  for i = 0, nslice - 1 do begin

     ;- velocity to extract, in pixel coordinates
     v = p.vcen - .95 * p.dv + 1.9 * p.dv * i / (nslice - 1)
     v_index = v

     ;- velocity in m/s
     vel = (v - sxpar(h, 'CRPIX3') + 1) * sxpar(h, 'CDELT3') + $
           sxpar(h, 'CRVAL3')

     slice = reform(co[*,*,v]) < 253B
     slice_raw = reform(co_raw[*,*,v])
     
     x = i mod ncol
     y = i / ncol
     x0 = xoff + x * dx
     y0 = yoff + (nrow - y - 1) * dy
     pos = [x0, y0 , $
            x0 + dx - xborder, $
            y0 + dy - yborder]
     indices, slice, x, y
     xyzadv, h, x, y, x * 0, a, d, junk
     a /= 15.

     ;- get formatted RA tick marks for the plot
     ra_ticks, a, tnum, tname, tval, number = 4

     tvimage, slice, /keep, pos = pos, /noi

     ;- bottom left plot gets axis information
     if i eq 6 then begin
        contour, slice_raw, a, d, pos = pos, /noerase, $
                 xsty = 1, ysty = 1, xra = [max(a), min(a)], $
                 xtit = textoidl('\alpha (J2000)'), $
                 ytit = textoidl('\delta (J2000)'), $
                 xtickname = textoidl(tname), xtickv = tval, $
                 xticks = tnum, xminor = 4, $
                 lev = lev
     endif else begin
        contour, slice_raw, a, d, pos = pos, /noerase, $
                 xsty = 1, ysty = 1, xra = [max(a), min(a)], $
                 lev = lev, $
                 xtickname = replicate(' ', tnum+1), $
                 ytickname = replicate(' ', 10), $
                 xticks = tnum, xminor = 4, xtickv = tval
     endelse

     ;- overplot candidate driver
     if src_ra[number] gt 0 then $
        oplot, src_ra[number]+[0,0], src_dec[number]+[0,0], psym = 8, symsize = 2, $
               color = 254

     ;- overplot model shell
     ;- R(v) = R_shell * sqrt(1 - ((v - v_cen)/v_exp)^2)
     contour, slice, pos = pos, /noerase, $
              xsty = 5, ysty = 5, $
              color = 255, c_thick = 2.0, /nodata
     rad = p.r * sqrt(1 - ((v - p.vcen) / p.dv)^2)
     oplot, p.ra + rad * xcircle, p.dec + rad * ycircle, $
            thick = 2, color = 255     

     ;- label velocity channel
     x = pos[0] + [.62, .98, .98, .62, .62] * (pos[2] - pos[0])
     y = pos[1] + [.87, .87, .98, .98, .87] * (pos[3] - pos[1])
     xy = convert_coord(x, y, /normal,/to_data)
     polyfill, xy[0,*], xy[1,*], color = 253
     oplot, xy[0,*], xy[1,*], color = 0, thick = 2
     xyouts, pos[0] + .96 * (pos[2] - pos[0]), $
             pos[1] + .89 * (pos[3] - pos[1]), $
             string(vel/1d3, format='(f0.1, " km/s")'), /norm, $
             color = 0, align = 1

  endfor
  device, /close
  set_plot, 'x'
end

pro driver

  channel_fig, 0, /co13
  channel_fig, 1, /co13
  for i = 2, 11 do channel_fig, i

  ;spawn, "for x in `ls channel_??.eps`; do epstopdf $x; done"
  ;spawn, "pdftk `ls channel_??.pdf` cat output channel_all.pdf"
end
