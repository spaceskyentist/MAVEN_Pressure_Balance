pro pbb_quicklook_forpublic
; Created 11/4/2024 SRS
; loads database of pressure balance boundary and generates plots for
; SRS et al. manuscript "Piercing the Martian Veil..."


;loadcv, 74, rgb_table=count_table,/noqual,/reverse
loadcv, 74, rgb_table=count_table,/noqual,/reverse
count_table[*,0]=!color.white
loadcv, 79, rgb_table=norm_table,/noqual
loadcv, 72, rgb_table=rdylbl_table, /noqual ; blue red4

plot_path = '/PressureBalance_Mars_SRS/'
plot_name_add = '_btot.png' ;'_btang.png' ;'_btot.png' ;
plot_title_add = ' Pb=Btot' ;' Pb=Bperp' ; add this to plot titles ;' Pb=Btot' ;

filename_db='full_reorg_db_2015-1-1_to_2021-2-28.sav'
load_pbb_db_for_analysis,filename_db=filename_db, database=db

filename_db='db_for_transitions_only_2015-01-01_to_2021-02-28.sav'
LOAD_ONELOCTRANS_PBB_DB, database=db_transoneloc, filename_db=filename_db

;number of transitions on an orbit
nt = abs(reform(median(db[*,29,*],dimension=1)))
iowt = where(nt gt 0., nowt) ; orbits with transition
iowot = where(nt eq 0., nowot) ; orbit segment without transition

;goto, MSOlocs
;goto, halekas

; ----------------------------- -----------------------------
; Local Time/ Altitude
; ----------------------------- -----------------------------
db_alt = reform(db[*,23,*], n_elements(db[*,23,*]))
db_lt = reform(db[*,24,*], n_elements(db[*,24,*]))
idb_finite = where(finite(db_alt),/null)
trans_loc = reform(db[*,29,iowt], n_elements(db[*,29,iowt]))
itrans = where(finite(trans_loc),/null)
trans_alt = reform(db[*,23,iowt], n_elements(db[*,23,iowt]))
trans_alt = trans_alt[itrans]
trans_lt = reform(db[*,24,iowt], n_elements(db[*,24,iowt]))
trans_lt = trans_lt[itrans]

xrange=[0,24] ;set binsize for all:
yrange=[200,800]
ybinsize = 10.
xbinsize = 1.
sstrans = create_grid(trans_lt, trans_alt,  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
sstot = create_grid(db_lt[idb_finite],  db_alt[idb_finite],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnorm = float(sstrans.counts)/float(sstot.counts) ; norm to total positions covered (above)
ssnorm[where(~finite(ssnorm),/null)] = 0.
ssnorm = ssnorm*100.D
;ssnorm_cb = ssnorm >30.
;ssnorm_cb = ssnorm_cb <70.


pntmp = plot_path+'lt_pbb/mvn_ltmso_alt_trans_viable_orbits'+plot_name_add
ptittmp = 'LT/Alt of Transitions !C b/n 200-800 km Altitude,'+plot_title_add ; Pb=abs(br)'
plot_2d_histogram, sstrans.xbins, sstrans.ybins, sstrans.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,400], rgbt=count_table, $
  plot_title=ptittmp, cb_title='Counts', xtitle='LT MSO [hr]', ytitle='Alt [km]'  ,$
  save_path_fn=pntmp

pntmp =  plot_path+'lt_pbb/mvn_ltmso_alt_all_viable_orbits'+plot_name_add
ptittmp = 'LT/Alt Visited on All Viable Orbits !C b/n 200-800 km Altitude,'+plot_title_add
plot_2d_histogram, sstot.xbins, sstot.ybins, sstot.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,400], rgbt=count_table, $
  plot_title=ptittmp, cb_title='Counts', xtitle='LT MSO [hr]', ytitle='Alt [km]',$
  save_path_fn=pntmp

pntmp = plot_path+'lt_pbb/mvn_ltmso_alt_trans_viable_orbits_norm2tot'+plot_name_add
ptittmp = 'LT/Alt  of Transitions - Normalized !C b/n 200-800 km Altitude,'+plot_title_add
plot_2d_histogram, sstrans.xbins, sstrans.ybins, ssnorm, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,400], rgbt=norm_table, $
  plot_title=ptittmp, cb_title='Percent', xtitle='LT MSO [hr]', ytitle='Alt [km]',$
  save_path_fn=pntmp



; Chi square fit
; --------------------------------------
; ; fit params from Chi Square fit - chisquareopt_conicfit.pro
theta_deg = [0:180]
theta_rad = theta_deg/!radeg
eccfit = 0.02
Lfit = 1.12
xffit = -0.0024
rfit = Lfit/(1.+eccfit*cos(theta_rad))
x_fit = rfit*cos(theta_rad) + xffit
y_fit = rfit*sin(theta_rad)


; Local hour on mars
theta_deg = [0:180]
lh = theta_Deg*(12D/180D)-12D ; 0 = midnight, 12=noon
;lh -= 24D*double(floor(lh/24D))     ; wrap to 0-24 range
lhdusk = lh+24D ; -12-0
lhdawn = abs(lh) ; 0-12 hr
rfull = [(rfit-1.)*3393., (rfit-1.)*3393.]
lh = [lhdusk, lhdawn]
ii = sort(lh)
lh=lh[ii]
rfull=rfull[ii]

onetrans_alt = reform(db_transoneloc[*,23])
onetrans_lt = reform(db_transoneloc[*,24,*])
ssonetrans = create_grid(onetrans_lt, onetrans_alt,  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
pntmp = plot_path+'lt_pbb/mvn_ltmso_alt_trans_viable_orbits_oneloctrans'+plot_name_add
ptittmp = 'LT/Alt  of Transitions - Using One Location of Trans!C b/n 200-800 km Altitude,'+plot_title_add
plot_2d_histogram, ssonetrans.xbins, ssonetrans.ybins, ssonetrans.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,400], rgbt=count_table, $
  plot_title=ptittmp, cb_title='Counts', xtitle='LT MSO [hr]', ytitle='Alt [km]',$
  save_path_fn=pntmp
p = plot(lh, rfull, linestyle='-', symbol='.', xtitle='LT', ytitle='Alt MSO [km]',/overplot,thick=4)
pntmp=plot_path+'lt_pbb/mvn_ltmso_alt_trans_onepostrans_viable_orbits_withfit'+plot_name_add
p.save, pntmp, border=10,/transparent,resolution=600


; ----------------------------- -----------------------------
; Lat/Lon
; ----------------------------- -----------------------------
db_lat = reform(db[*,22,*], n_elements(db[*,22,*]))
db_lon = reform(db[*,21,*], n_elements(db[*,21,*]))
idb_finite = where(finite(db_alt),/null)
trans_loc = reform(db[*,29,iowt], n_elements(db[*,29,iowt]))
itrans = where(finite(trans_loc),/null)
trans_lat = reform(db[*,22,iowt], n_elements(db[*,22,iowt]))
trans_lat = trans_lat[itrans]
trans_lon = reform(db[*,21,iowt], n_elements(db[*,21,iowt]))
trans_lon = trans_lon[itrans]

;;; Get the crustal fields
restore, '/Users/sksh3793/maven/plot/BrainOrbitGeometryPlotFiles/br_360x180_pc.sav'
mapcolor = colorscale( br, mindat=-70, maxdat=70, mincol=1, maxcol=254 )
tmp = where(br eq 0, tmpcnt)
if tmpcnt ne 0 then mapcolor[tmp] = 255
br[where(abs(br) lt 10.)] = 0. ; only plot crustal fields with magnitude greater than 10 nT

xrange=[0,360.]
yrange=[-90.,90.]
binsize= 10. ; degrees
ss = create_grid(trans_lon, trans_lat, miny=min(yrange),maxy=max(yrange), ybinsize=binsize, $ ; Trans locs
  minx=min(xrange),maxx=max(xrange), xbsize=binsize)

pntmp = plot_path+'latlon_pbb/mvn_lat_lon_trans_viable_orbits_pbdomhigh'+plot_name_add
ptittmp = 'Lat/Lon Locations Of Transitions B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
plot_2d_histogram, ss.xbins, ss.ybins, ss.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[900,600], rgbt=count_table, $
  plot_title=Ptittmp, cb_title='Counts', xtitle='Longitude [deg]', ytitle='Latitude [deg]', $
  /plotmars, change_aspect_ratio=0, $
  save_path_fn=pntmp
lattmp = [-90:90-1]
lontmp = [0:360-1]
c_cf = contour(br,lontmp, lattmp,/overplot, rgb_table=70, color='k')  ;, c_value=v) ;color='k')
c_cf.save, plot_path+'latlon_pbb/mvn_lat_lon_trans_viable_orbits_pbdomhigh'+plot_name_add


; get index of periapsis on orbits without transition
dbtmp = db[*,*,iowot] ; without transitions
ntmp = n_elements(iowot) ; number of orbits
sztmp = size(dbtmp,/dimensions)
alttmp = reform(dbtmp[*,23,*])
imin = fltarr(ntmp)*!values.f_nan ; find periapsis index
for ii=0,ntmp-1 do imin[ii] = min_index(alttmp[*,ii])
peri_lat =[]
peri_lon = []
for ii=0,ntmp-1 do peri_lat = [peri_lat, dbtmp[imin[ii],22,ii]]
for ii=0,ntmp-1 do peri_lon = [peri_lon, dbtmp[imin[ii],21,ii]]

ptittmp = 'Lat/Lon Locations of MAVEN at 200 km !C Orbits without Transitions B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
ss = create_grid(peri_lon, peri_lat, miny=min(yrange),maxy=max(yrange), ybinsize=binsize, $ ; MSO locations visited on orbits with transitions
  minx=min(xrange),maxx=max(xrange), xbsize=binsize)
pntmp =plot_path+'latlon_pbb/mvn_lat_lon_NOtrans_viable_orbits_pbdomhigh'+plot_name_add
plot_2d_histogram, ss.xbins, ss.ybins, ss.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[900,600], rgbt=count_table, $
  plot_title=Ptittmp, cb_title='Counts', xtitle='Longitude [deg]', ytitle='Latitude [deg]', $
  /plotmars, change_aspect_ratio=0, $
  save_path_fn=pntmp
lattmp = [-90:90-1]
lontmp = [0:360-1]
c_cf = contour(br,lontmp, lattmp,/overplot, rgb_table=70, color='k')  ;, c_value=v) ;color='k')
c_cf.save, plot_path+'latlon_pbb/mvn_lat_lon_NOtrans_viable_orbits_pbdomhigh'+plot_name_add


; averaged trans loc
; LAT/LON
lon = reform(db_transoneloc[*,21])
lat = reform(db_transoneloc[*,22])
ptittmp = 'Lat/Lon Locations Of Transitions B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
ss = create_grid(lon, lat, miny=min(yrange),maxy=max(yrange), ybinsize=binsize, $ ; MSO locations visited on orbits with transitions
  minx=min(xrange),maxx=max(xrange), xbsize=binsize)
pntmp =  plot_path+'latlon_pbb/mvn_lat_lon_trans_viable_orbits_pbdomhigh_averagedtranslocation'+plot_name_add
plot_2d_histogram, ss.xbins, ss.ybins, ss.counts, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[900,600], rgbt=count_table, $
  plot_title=Ptittmp, cb_title='Counts', xtitle='Longitude [deg]', ytitle='Latitude [deg]', $
  /plotmars, change_aspect_ratio=0, $
  save_path_fn=pntmp
lattmp = [-90:90-1]
lontmp = [0:360-1]
c_cf = contour(br,lontmp, lattmp,/overplot, rgb_table=70, color='k')  ;, c_value=v) ;color='k')
c_cf.save, plot_path+'latlon_pbb/mvn_lat_lon_trans_viable_orbits_pbdomhigh_averagedtranslocation'+plot_name_add




; ----------------------------- -----------------------------
; MSO Locations
; ----------------------------- -----------------------------
MSOlocs:

db_Xmso = reform(db[*,16,*], n_elements(db[*,16,*]))
db_Ymso = reform(db[*,17,*], n_elements(db[*,17,*]))
db_Zmso = reform(db[*,18,*], n_elements(db[*,18,*]))
idb_finite = where(finite(db_Xmso),/null)
trans_loc = reform(db[*,29,iowt], n_elements(db[*,29,iowt]))
itrans = where(finite(trans_loc),/null)
trans_Xmso = reform(db[*,16,iowt], n_elements(db[*,16,iowt]))
trans_Xmso = trans_Xmso[itrans]
trans_Ymso = reform(db[*,17,iowt], n_elements(db[*,17,iowt]))
trans_Ymso = trans_Ymso[itrans]
trans_Zmso = reform(db[*,18,iowt], n_elements(db[*,18,iowt]))
trans_Zmso = trans_Zmso[itrans]

xrange=[-1.5,1.5] ;set binsize for all:
yrange=xrange
zrange=zrange
ybinsize = 0.1 ;0.05 ;0.03
xbinsize =  0.1 ;0.05 ;0.03
upper_cb_cut = 6.

;XY
sstransXY = create_grid(trans_Xmso, trans_Ymso,  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
sstotXY = create_grid(db_Xmso[idb_finite],  db_Ymso[idb_finite],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnormXY = float(sstransXY.counts)/float(sstotXY.counts) ; norm to total positions covered (above)
ssnormXY[where(~finite(ssnormXY),/null)] = 0.
ssnormXY[where(ssnormXY eq 0.,/null)] = !values.f_nan
ssnormXY[where(ssnormXY ge 0.2,/null)] = !values.f_nan
ssnormXY_cb = ssnormXY*100.D
;ssnormXY_cb = ssnormXY >2.
ssnormXY_cb = ssnormXY_cb < upper_cb_cut

;XY positive Z
sstransXYposZ = create_grid(trans_Xmso[where(trans_Zmso ge 0.,/null)], trans_Ymso[where(trans_Zmso ge 0.,/null)],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
xtmp = db_Xmso[idb_finite]
ytmp = db_Ymso[idb_finite]
ztmp = db_Zmso[idb_finite]
sstotXYposZ = create_grid(xtmp[where(ztmp ge 0.,/null)],  ytmp[where(ztmp ge 0.,/null)],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnormXYposZ = float(sstransXYposZ.counts)/float(sstotXYposZ.counts) ; norm to total positions covered (above)
ssnormXYposZ[where(~finite(ssnormXYposZ),/null)] = 0.
ssnormXYposZ[where(ssnormXYposZ eq 0.,/null)] = !values.f_nan
ssnormXYposZ[where(ssnormXYposZ ge 0.2,/null)] = !values.f_nan
ssnormXYposZ_cb = ssnormXYposZ*100.D
;ssnormXY_cb = ssnormXY >2.
ssnormXYposZ_cb = ssnormXYposZ_cb < upper_cb_cut

;XY negative Z
sstransXYnegZ = create_grid(trans_Xmso[where(trans_Zmso lt 0.,/null)], trans_Ymso[where(trans_Zmso lt 0.,/null)],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
xtmp = db_Xmso[idb_finite]
ytmp = db_Ymso[idb_finite]
ztmp = db_Zmso[idb_finite]
sstotXYnegZ = create_grid(xtmp[where(ztmp lt 0.,/null)],  ytmp[where(ztmp lt 0.,/null)],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnormXYnegZ = float(sstransXYnegZ.counts)/float(sstotXYnegZ.counts) ; norm to total positions covered (above)
ssnormXYnegZ[where(~finite(ssnormXYnegZ),/null)] = 0.
ssnormXYnegZ[where(ssnormXYnegZ eq 0.,/null)] = !values.f_nan
ssnormXYnegZ[where(ssnormXYnegZ ge 0.2,/null)] = !values.f_nan
ssnormXYnegZ_cb = ssnormXYnegZ*100.D
;ssnormXY_cb = ssnormXY >2.
ssnormXYnegZ_cb = ssnormXYnegZ_cb < upper_cb_cut

;XZ
sstransXZ = create_grid(trans_Xmso, trans_Zmso,  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
sstotXZ = create_grid(db_Xmso[idb_finite],  db_Zmso[idb_finite],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnormXZ = float(sstransXZ.counts)/float(sstotXZ.counts) ; norm to total positions covered (above)
ssnormXZ[where(~finite(ssnormXZ),/null)] = 0.
ssnormXZ[where(ssnormXZ eq 0.,/null)] = !values.f_nan
ssnormXZ_cb = ssnormXZ*100.D
;ssnormXZ_cb = ssnormXZ >2.
ssnormXZ_cb = ssnormXZ_cb < upper_cb_cut

;YZ
sstransYZ = create_grid(trans_Ymso, trans_Zmso,  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
sstotYZ = create_grid(db_Ymso[idb_finite],  db_Zmso[idb_finite],  miny=min(yrange),maxy=max(yrange), ybinsize=ybinsize, $ ;
  minx=min(xrange),maxx=max(xrange), xbsize=xbinsize)
ssnormYZ = float(sstransYZ.counts)/float(sstotYZ.counts) ; norm to total positions covered (above)
ssnormYZ[where(~finite(ssnormYZ),/null)] = 0.
ssnormYZ[where(ssnormYZ eq 0.,/null)] = !values.f_nan
ssnormYZ_cb = ssnormYZ*100.D
;ssnormYZ_cb = ssnormYZ >2.
ssnormYZ_cb = ssnormYZ_cb < upper_cb_cut

loadcv, 64, rgb_table=count_table_mso,/noqual ,/reverse


; PLOTS!
;---------------------

;XY
ptittmp = 'XY MSO Locations of Transitions - Normalized !C B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
pntmp =plot_path+'mso_pbb/mvn_mso_trans_viable_orbits_pbdom_norm2tot_x_y'+plot_name_add
plot_2d_histogram, sstransxy.xbins, sstransxy.ybins, ssnormXY_cb, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,800], rgbt=norm_table, $ ;rgbt=count_table_mso, $
  plot_title=ptittmp, cb_title='Percent', xtitle='X MSO', ytitle='Y MSO', $
  /plotmars,save_path_fn=pntmp, extraprops={buffer:1}
;c = contour(sstransxy.counts, sstransxy.xbins, sstransxy.ybins,/overplot)
;c = contour(ssnormXY*100., sstransxy.xbins, sstransxy.ybins,/overplot,rgb_table=0, n_levels=2) ;,/fill)

;XY positive Z
ptittmp = 'XY MSO Locations of Transitions in +Z MSO hemisphere - Normalized !C B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
pntmp =plot_path+'mso_pbb/mvn_mso_trans_viable_orbits_pbdom_norm2tot_x_y_posz'+plot_name_add
plot_2d_histogram, sstransxyposZ.xbins, sstransxyposZ.ybins, ssnormXYposZ_cb, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,800], rgbt=norm_table, $
  plot_title=ptittmp, cb_title='Percent', xtitle='X MSO', ytitle='Y MSO', $
  /plotmars;,save_path_fn=pntmp, extraprops={buffer:1}

;XY negative Z
ptittmp = 'XY MSO Locations of Transitions in -Z MSO hemisphere - Normalized !C B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
pntmp =plot_path+'mso_pbb/mvn_mso_trans_viable_orbits_pbdom_norm2tot_x_y_negz'+plot_name_add
plot_2d_histogram, sstransxynegZ.xbins, sstransxynegZ.ybins, ssnormXYnegZ_cb, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,800], rgbt=norm_table, $
  plot_title=ptittmp, cb_title='Percent', xtitle='X MSO', ytitle='Y MSO', $
  /plotmars;,save_path_fn=pntmp, extraprops={buffer:1}

;XZ
ptittmp = 'XZ MSO Locations of Transitions - Normalized !C B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
pntmp =plot_path+'mso_pbb/mvn_mso_trans_viable_orbits_pbdom_norm2tot_x_z'+plot_name_add
plot_2d_histogram, sstransxz.xbins, sstransxz.ybins, ssnormXZ_cb, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,800], rgbt=norm_table, $ ;rgbt=count_table_mso, $
  plot_title=ptittmp, cb_title='Percent', xtitle='X MSO', ytitle='Z MSO', $
  /plotmars,save_path_fn=pntmp, extraprops={buffer:1}

;YZ
ptittmp = 'YZ MSO Locations of Transitions - Normalized !C B/N 200-800 km Altitude, Pb dominant at high alts,'+plot_title_add
pntmp =plot_path+'mso_pbb/mvn_mso_trans_viable_orbits_pbdom_norm2tot_y_z'+plot_name_add
plot_2d_histogram, sstransyz.xbins, sstransyz.ybins, ssnormYZ_cb, $
  xrange=xrange, yrange=yrange, binsize=binsize, winsize=[800,800], rgbt=norm_table, $ ;rgbt=count_table_mso, $
  plot_title=ptittmp, cb_title='Percent', xtitle='Y MSO', ytitle='Z MSO', $
  /plotmars,save_path_fn=pntmp, extraprops={buffer:1}




end