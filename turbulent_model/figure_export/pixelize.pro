pro pixelize, p, h
  advxyz, h, p.ra, p.dec, p.vcen, x, y, z
  p.r /= abs(sxpar(h, 'cdelt1'))
  p.dr /= abs(sxpar(h, 'cdelt1'))
  p.dv /= abs(sxpar(h, 'cdelt3'))
  p.ra = x
  p.dec = y
  p.vcen = z
end

  
  
