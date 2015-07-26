#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Communicate w/ a Decathlon Keymaze 500/700 devices
#-----------------------------------------------------------------------------
# @author Emmanuel Blot <manu.blot@gmail.com> (c) 2009
# @license MIT License, see LICENSE file
#-----------------------------------------------------------------------------

from __future__ import with_statement
from optparse import OptionParser
from db import KeymazeCache
from keymaze import KeymazePort
import datetime
import logging
import math
import os
import re
import time
import sys
import glob
import serial

def show_trackpoints_catalog(cat):
    sorted_cat = list(cat)
    sorted_cat.sort(key=lambda e: e['id'])
    print " #   Day         Start  End    Duration  Distance  AltMin   AltMax    Laps"
    for tpent in sorted_cat:
        stime = time.localtime(tpent['start'])
        print ' %02d  %s  %s  %s    %s  %6.2fkm %s  %s  %s' % \
            (tpent['id'],
             time.strftime('%Y-%m-%d', stime),
             time.strftime('%H:%M', stime),
             time.strftime('%H:%M', 
                           time.localtime(tpent['start']+tpent['time'])),
             tpent['duration'] and \
                time.strftime('%Hh%Mm', 
                              time.localtime(tpent['duration']-3600)) or
                '     -',
             tpent['distance']/1000.0,
             tpent['altmin'] and '%6dm' % tpent['altmin'] or '      -', 
             tpent['altmax'] and '%6dm' % tpent['altmax'] or '      -',
             tpent['laps'] and '%6d' % tpent['laps'] or '      -')

def parse_trim(trim_times):
    tcre = re.compile(r'^(?P<r>[+-])?'
                      r'(?:(?P<h>\d\d):(?=\d\d:))?'
                      r'(?:(?P<m>\d\d):)?'
                      r'(?P<s>\d\d)$')
    values = []
    for trim in trim_times:
        mo = tcre.match(trim)
        if not mo:
            raise AssertionError('Invalid trim format "%s"' % trim)
        seconds = 60*int(mo.group('h') or 0)
        seconds = 60*(seconds + int(mo.group('m') or 0))
        seconds += int(mo.group('s') or 0)
        values.append((mo.group('r'), seconds))
    return values

def trim_trackpoints(track, tp, trims):
    if not trims:
        return tp
    start = track['start']
    end = start + track['time']
    if len(trims) > 1:
        tstart = trims[0]
        tend = trims[1]
    else:
        tstart = trims[0]
        tend = ('', 0)
    if tstart[0] == '+':
        t_start = start+tstart[1]
    elif tstart[0] == '-':
        t_start = end-tstart[1]
    else:
        # not yet implemented
        t_start = start
    if tend[0] == '+':
        t_end = start+tend[1]
    elif tend[0] == '-':
        t_end = end-tend[1]
    else:
        # not yet implemented
        t_end = end
    ttp = []
    pt = start*10
    t_start *= 10
    t_end *= 10
    for p in tp:
        pt += p[-1] # delta
        if t_start <= pt <= t_end:
            ttp.append(p)
    return ttp
    
def haversine(pt1, pt2):
    (lat1, lon1) = map(math.radians, pt1[0:2])
    (lat2, lon2) = map(math.radians, pt2[0:2])
    R = 6371.0*1000 # m
    # Mean radius        6,371.0 km
    # Equatorial radius  6,378.1 km
    # Polar radius       6,356.8 km
    dlat = lat2-lat1
    dlon = lon2-lon1 
    a = math.sin(dlat/2) * math.sin(dlat/2) + \
        math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2) * math.sin(dlon/2) 
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) 
    d = R * c
    return d
    
def cartesian(point):
    lat = math.radians(point[0])
    lon = math.radians(point[1])
    R = 6371.0*1000 # m
    h = R + point[2]
    x = h * math.cos(lat) * math.cos(lon)
    y = h * math.cos(lat) * math.sin(lon)
    z = h * math.sin(lat)
    return (x,y,z)

def dotproduct(a,b):
    return sum((a[0]*b[0],a[1]*b[1],a[2]*b[2]))
    
def delta(c1, c2, c3):
    (x1,y1,z1) = c1
    (x2,y2,z2) = c2
    (x3,y3,z3) = c3
    u = ((x2-x1),(y2-y1),(z2-z1))
    v = ((x3-x2),(y3-y2),(z3-z2))
    uv = dotproduct(u,v)
    lu = math.sqrt(dotproduct(u,u))
    lv = math.sqrt(dotproduct(v,v))
    try:
        theta = math.degrees(math.acos(uv/(lu*lv)))
    except ZeroDivisionError:
        theta = 0
    return (lu, theta, lv)

def optimize(points, angle=0):
    tpoints = [(tp[0]/1000000.0, tp[1]/1000000.0, tp[2], tp[3], tp[4], tp[5]) \
        for tp in points]
    if angle == 0:
        return tpoints
    queue = [tpoints[0], tpoints[0]]
    opt = list(queue)
    for tp in tpoints:
        # d = haversine(last[1], tp)
        da = delta(cartesian(queue[0]), cartesian(queue[1]), cartesian(tp))
        queue.pop(0)
        queue.append(tp)
        #print da
        if da[1] > angle:
            opt.append(tp)
    return opt

def get_available_serial_ports():
    # Code from Thomas @ Stackoverflow
    # http://stackoverflow.com/questions/12090503/listing-available-com-ports-with-python
    """Lists serial ports

    :raises EnvironmentError:
        On unsupported or unknown platforms
    :returns:
        A list of available serial ports
    """
    if sys.platform.startswith('win'):
        ports = ['COM' + str(i + 1) for i in range(256)]

    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this is to exclude your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')

    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty.*')

    else:
        raise EnvironmentError('Unsupported platform')

    result = []
    for port in ports:
        try:
            s = serial.Serial(port)
            s.close()
            result.append(port)
        except (OSError, serial.SerialException):
            pass
    return result


class Keymaze(object):
    def __init__(self, dbpath = None, port=None, offline_mode=False):
        log = logging.getLogger('pykmaze')
        log.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter("%(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        log.addHandler(ch)
        self.log = log
        
        if port is None:
            serial_ports = get_available_serial_ports()
            if len(serial_ports) == 0:
                print 'No serial port detected, going to offline mode'
                offline_mode = True
            else:
                print 'Available serial port:'
                print serial_ports
                port = serial_ports[0]
                print 'Taking %s' % port
        
        if dbpath is None:
            dbpath = self._get_default_dbpath()
        self.dbpath = dbpath
        
        if not offline_mode:
            self.set_online_mode(port)
        else:
            self.set_offline_mode()
    
    def set_offline_mode(self):
        self.init_keymaze(None)
        self.current_mode = 'offline'
    
    def set_online_mode(self, port):
        keymaze_port = KeymazePort(self.log, port)
        self.init_keymaze(keymaze_port)
        self.current_mode = 'online'
        
    def init_keymaze(self, keymaze_port):
        self.cache = KeymazeCache(self.log, self.dbpath, keymaze_port)
        info = self.cache.get_information()
        self.device = self.cache.get_device(info['serialnumber'])
        self.tpcat = self.cache.get_trackpoint_catalog(self.device)
    
    def _get_default_dbpath(self):
        dbname = 'pykmaze.sqlite'
        dbpath = 'HOME' in os.environ and os.environ['HOME'] or '.'
        if sys.platform.lower() in ('darwin'):
            dbpath = os.path.join(dbpath, 'Library', 'Application Support', 
                                  'Pykmaze', dbname)
        else:
            dbpath = os.path.join(dbpath, '.pykmaze', dbname)
        return dbpath
    
    def print_info(self):
        print ' Device: %s' % self.info['name']
        print ' Owner:  %s' % self.info['user']
        print ' S/N:    %s' % self.info['serialnumber']
        print ''
    
    def sync(self):
        if self.current_mode == 'offline':
            raise AssertionError('Cannot sync from device in offline mode')
        reload_cache = False
        for tp in self.tpcat:
            if None in (tp['altmin'], tp['altmax']):
                self.log.info('Should load sync track %d from device' % \
                         (int(tp['id'])))
                self.cache.get_trackpoints(self.device, tp)
                reload_cache = True
        if reload_cache:
            self.tpcat = self.cache.get_trackpoint_catalog(self.device)
    
    def show_catalog(self):
        show_trackpoints_catalog(self.tpcat)
        print ''
    
    def get_all_tracks(self):
        tracks = [tp['track'] for tp in self.tpcat]
        self.log.debug('Tracks %s' % tracks)
        return tracks
        
    def get_trackinfo_from_name(self, trackname):
        trackname = int(trackname)
        for tp in self.tpcat:
            if int(tp['id']) == trackname:
                return tp
        raise AssertionError('Track "%s" does not exist' % \
                                         trackname)
    
    def get_trackpoints(self, trackinfo):
        self.log.info('Recovering trackpoints for track %u' % trackinfo['id'])
        tpoints = self.cache.get_trackpoints(self.device, trackinfo)
        return tpoints
    
    def get_trackinfo(self, track):
        return filter(lambda x: x['id'] == track['id'], 
                                        self.tpcat)[0]
    
    def trim_trackpoints(self, track_info, tpoints, trim):
        trims = parse_trim(trim.split(','))
        self.log.info('All points: %d' % len(tpoints))
        tpoints = trim_trackpoints(track_info, tpoints, trims)
        self.log.info('Filtered points: %d' % len(tpoints))
        return tpoints
    
    def export_track(self, track, filename, filetype, trim=None,
                     prepend_datetime=False, zoffset=0, mode='default'):
        # filetype 'gpx', 'kmz' or 'km'
        if filetype not in ['kmz', 'kml', 'gpx']:
            raise AssertionError('Unsupported filetype %s' % filetype)
        
        track_info = self.get_trackinfo(track)
        tpoints = self.get_trackpoints(track_info)
        if prepend_datetime:
            stime = time.localtime(track_info['start'])
            filename = os.path.join(os.path.dirname(filename),
                       time.strftime('%Y-%m-%d-%H-%M-%S-', stime) + \
                       os.path.basename(filename))
        
        is_km = filetype in ['kmz', 'kml']
        if is_km:
            from kml import KmlDoc
        if filetype == 'gpx':
            from gpx import GpxDoc

        if trim:
            tpoints = trim_trackpoints(track_info, tpoints, trim)

        optpoints = optimize(tpoints, 0)
        self.log.info('Count: %u, opt: %u', 
                  len(tpoints), len(optpoints))
        if is_km:
            kml = KmlDoc(os.path.splitext(os.path.basename(filename))[0])
            kml.add_trackpoints(optpoints, int(zoffset), 
                                extrude='air' not in mode)
        if filetype == 'kmz':
            import zipfile
            import cStringIO as StringIO
            out = StringIO.StringIO()
            kml.write(out)
            out.write('\n')
            z = zipfile.ZipFile(filename, 'w', 
                                zipfile.ZIP_DEFLATED)
            z.writestr('doc.kml', out.getvalue())
        if filetype == 'kml':
            with open(filename, 'wt') as out:
                kml.write(out)
                out.write('\n')
        
        if filetype == 'gpx':
            gpx = GpxDoc(os.path.splitext( \
                         os.path.basename(filename))[0],
                         track_info['start'])
            gpx.add_trackpoints(optpoints, int(zoffset))
            with open(filename, 'wt') as out:
                gpx.write(out)
                out.write('\n')
    
    def export_all_tracks(self, basename, filetype):
        for track in self.get_all_tracks():
            self.export_track(track, basename, filetype, prepend_datetime=True)

if __name__ == '__main__':
    dbname = 'pykmaze.sqlite'
    dbpath = 'HOME' in os.environ and os.environ['HOME'] or '.'
    if sys.platform.lower() in ('darwin'):
        dbpath = os.path.join(dbpath, 'Library', 'Application Support', 
                              'Pykmaze', dbname)
    else:
        dbpath = os.path.join(dbpath, '.pykmaze', dbname)
    modes = ('default', 'air')
    usage = 'Usage: %prog [options]\n' \
            '   Keymaze 500-700 communicator'
    optparser = OptionParser(usage=usage)
    optparser.add_option('-p', '--port', dest='port',
                         default=None,
                         help='Serial port name')
    optparser.add_option('-k', '--kml', dest='kml',
                         help='Export to KML, output file name')
    optparser.add_option('-K', '--kmz', dest='kmz',
                         help='Export to KMZ, output file name')
    optparser.add_option('-x', '--gpx', dest='gpx',
                         help='Export to GPX, output file name')
    optparser.add_option('-T', '--trim', dest='trim',
                         help='Trim a track start[,end] with [+-]hh:mm:ss')
    optparser.add_option('-z', '--zoffset', dest='zoffset', default='0',
                         help='Offset to add on z-axis (meters)')
    optparser.add_option('-s', '--storage', dest='storage', 
                         default=dbpath,
                         help='Specify path for data storage (default: %s)' \
                                % dbpath)
    optparser.add_option('-o', '--offline', dest='offline', 
                         action='store_true',
                         help='Offline (used cached information)')
    optparser.add_option('-S', '--sync', dest='sync', 
                         action='store_true',
                         help='Load all new tracks from device')
    optparser.add_option('-f', '--force', dest='force', 
                         action='store_true',
                         help='Force reload from device')
    optparser.add_option('-i', '--info', dest='info', 
                         action='store_true',
                         help='Show owner information')
    optparser.add_option('-c', '--catalog', dest='catalog', 
                         action='store_true',
                         help='Show track catalog')
    optparser.add_option('-t', '--track', dest='track', 
                         help='Retrieve trackpoint for specified track')
    optparser.add_option('-m', '--mode', dest='mode', choices=modes,
                         help='Use show mode among [%s]' % ','.join(modes),
                         default=modes[0])
    
    (options, args) = optparser.parse_args(sys.argv[1:])
    
    
    try:
        if options.force and options.offline:
            raise AssertionError('Force and offline modes are mutually '
                                 'exclusive')
        keymaze = Keymaze(dbpath=options.storage, port=options.port, offline_mode=options.offline)

        if options.info:
            keymaze.print_info()
        
        if options.sync:
            keymaze.sync()
                    
        if options.catalog:
            keymaze.show_catalog()
        
        if not options.track or options.track == 'all':
            tracks = keymaze.tpcat
            prepend_datetime = True
        else:
            tracks = [keymaze.get_track_from_name(options.track)]
            prepend_datetime = False

        for track in tracks:
            if options.kml:
                keymaze.export_track(track, options.kml, 'kml',
                                     options.trim, prepend_datetime,
                                     options.zoffset, options.mode)
            if options.kmz:
                keymaze.export_track(track, options.kmz, 'kmz',
                                     options.trim, prepend_datetime,
                                     options.zoffset, options.mode)
            if options.gpx:
                keymaze.export_track(track, options.gpx, 'gpx',
                                     options.trim, prepend_datetime,
                                     options.zoffset, options.mode)
                        
    except AssertionError, e:
        print >> sys.stderr, 'Error: %s' % e[0]
