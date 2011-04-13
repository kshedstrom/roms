load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  cdf_file = addfile("chukchi_sta.nc", "r")

  temp1 = cdf_file->temp(:,4,:)
  t0 = transpose(temp1)
  salt1 = cdf_file->salt(:,4,:)
  s0 = transpose(salt1)
  u1 = cdf_file->salt(:,4,:)
  u0 = transpose(u1)

  sc_r = cdf_file->s_rho(:)
  Cs_r = cdf_file->Cs_r(:)
  hc = cdf_file->hc
  h = 10.0

  z = hc*sc_r + (h-hc)*Cs_r

  time = cdf_file->ocean_time(:)
  time = (time/86400) - 31052

;  wks = gsn_open_wks("x11","sta3")  ; Open an X11 workstation.
  wks = gsn_open_wks("ncgm","sta3") ; Open an NCGM workstation.

;---------- Begin plots -----------------------------------------

  resources                    = True          ; Indicate you want to
                                               ; set some resources.

  resources@xyLineColors  = (/10,12,16,4,14,6,3,17,9,15,11,5,13,7,2,\
                              10,12,16,4,14,6,3,17,9,15,11,5,13,7,2/)
  resources@xyLineThicknesses   = 3

;  resources@tiMainString = "Susitna 1"
  plot = gsn_xy(wks,time,t0,resources)
  plot = gsn_xy(wks,time,s0,resources)
  plot = gsn_xy(wks,time,u0,resources)
;  resources@tiMainString = "Susitna 2"
;  plot = gsn_xy(wks,time,t3,resources)
;  plot = gsn_xy(wks,time,s3,resources)

;  resources@tiMainString = "Stikine 1"
;  plot = gsn_xy(wks,time,t4,resources)
;  plot = gsn_xy(wks,time,s4,resources)
;  resources@tiMainString = "Stikine 2"
;  plot = gsn_xy(wks,time,t9,resources)
;  plot = gsn_xy(wks,time,s9,resources)

;  resources@sfYCStartV = 0
;  resources@sfYCEndV = -10
;  resources@sfYArray = z
;  resources@cnFillOn = True
;  resources@tiYAxisString = "Distance from Surface"
;  resources@tiXAxisString = "Time since Initialization (days)"
;  resources@cnFillColors  = (/10,12,16,4,14,6,3,17,9,15,11,5,13,7,2,\
;                              10,12,16,4,14,6,3,17,9,15,11,5,13,7,2/)
;
;  plot = gsn_contour(wks, t0, resources)
;  plot = gsn_contour(wks, s0, resources)
;
;  delete(plot)       ; Clean up.
;  delete(resources)
end
