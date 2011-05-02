load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  cdf_file = addfile("chukchi_sta.nc", "r")

;  temp1 = cdf_file->temp(2700:3000,180:187,41)
  temp1 = cdf_file->temp(:,:,1)
  t1 = transpose(temp1)
;  salt1 = cdf_file->salt(2700:3000,180:187,41)
  salt1 = cdf_file->salt(:,:,1)
  s1 = transpose(salt1)
;  uv1 = cdf_file->u(2700:3000,180:187,41)
  uv1 = cdf_file->u(:,:,1)
  u1 = transpose(uv1)
;  time = cdf_file->ocean_time(2700:3000)
  time = cdf_file->ocean_time
  time = (time/86400) - 31052.

  wks = gsn_open_wks("x11","gsun01n")  ; Open an X11 workstation.

;---------- Begin third plot -----------------------------------------

  resources                    = True          ; Indicate you want to
                                               ; set some resources.

;  resources@xyLineColors        = (/8,9,10,11,12,13,14,14,16/)  ; Define line colors.
  resources@xyLineThicknesses   = 2
                                               ; (1.0 is the default).

  resources@tiMainString = "Bristol Bay"
  plot = gsn_xy(wks,time,t1,resources)
  plot = gsn_xy(wks,time,s1,resources)
  plot = gsn_xy(wks,time,u1,resources)

  delete(plot)       ; Clean up.
  delete(resources)
end
