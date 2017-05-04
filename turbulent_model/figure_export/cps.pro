function cps, number

  ra = [ten(3, 27, 35), ten(3, 27, 39), ten(3, 29, 25), $
        ten(3, 35, 25), ten(3, 40, 07), ten(3, 41, 24), $
        ten(3, 42, 17), ten(3, 44, 10), ten(3, 44, 20), $
        ten(3, 44, 34), ten(3, 44, 52), ten(3, 47, 46)] * 15.
  dec = [ten(31, 17, 00), ten(31, 04, 30), ten(31, 25, 10), $
         ten(31, 09, 00), ten(31, 50, 20), ten(31, 53, 40), $
         ten(32, 06, 40), ten(32, 17, 00), ten(32, 08, 30), $
         ten(32, 10, 00), ten(32, 18, 30), ten(33, 03, 40)]

  ;- radius, dr in degrees
  radius = [7.2, 5.8, 2.0, 13., 37, 7.8, 30, 13, 2.6, 2, 2.2, 6] / 60.
  dr = [2.5, 3.7, 1.3, 10, 4, 2.7, 5, 2, 1.4, 1, 1.3, 4] / 60.

  ;- vel, vcen in m/s
  ;- vel is cap-to-midplane
  vel = [2.5, 2.5, 1, 5, 5, 3, 3, 2.5, 2, 2, 2, 2.5] * 1000
;  vcen = [6.97, 5.84, 8.04, 4.39, 8.84, 8.17, 9.36, $
;          9.04, 8.90, 7.7, 8.3, 9.9] * 1000
  vcen = [7.4, 7.2, 7.85, 1.99, $
          5.04, 7.5, 6.2, 9.2, 9.1, $
          7.65, 7.34, 9.3] * 1000
  rec = {ra:0., dec:0., r:0., dr:0., vcen:0., dv:0.}
  
  result = replicate(rec, 12)
  result.ra = ra
  result.dec = dec
  result.r = radius
  result.dr = dr
  result.vcen = vcen
  result.dv = vel

  if n_elements(number) ne 0 then return, result[number]
  return, result
end
