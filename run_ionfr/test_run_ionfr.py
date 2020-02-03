from run_ionfr import run_ionfr

ra = '08h37m05.6s'
dec = '06d10m14.5s'
date = '2011-10-20'
lat = '52d54m54.6sn'
lon = '6d52m11.7se'
source = 'xx'

run_ionfr(ra, dec, source, date, sign= '+', convertradec = False, lat=lat, lon=lon)