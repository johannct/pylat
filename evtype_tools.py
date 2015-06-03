from IPython.parallel import require,Client
import pyfits
from IPython.display import clear_output
import os, subprocess

ALL_TYPES=["FRONT","BACK"]+["%s%d"%("PSF",i) for i in range(4)]+["%s%d"%("EDISP",i) for i in range(4)]
TYPES_DICT={"FB":[1,2],"PSF":[4,8,16,32],"EDISP":[64,128,256,512]}
P8R2_CLASS_DICT = {"TRANSIENT020E":8,"TRANSIENT020":16,"TRANSIENT010E":32,"TRANSIENT0101":64,"SOURCE":128,"CLEAN":256,"ULTRACLEAN":512,"ULTRACLEANVETO":1024,"TRANSIENT100S":32768,"TRANSIENT015S":65536}


def build_xml(workdir, xml_basefile, evtype):
    output = subprocess.check_output(['grep','txt',xml_basefile])
    isoname = output[output.find('iso'):].split('"')[0]

    new_file = os.path.join(workdir,xml_basefile.replace('.xml','_%s.xml'%evtype))
    output = os.system('cp 3Cs_3FGL_v14_15deg_withdiff.xml %s'%new_file)
    if not output==0:
        raise "problem with copying %s"%xml_basefile

    output = os.system("sed -i 's/%s/iso_P8R2_SOURCE_V6_%s_v06.txt/g' %s"%(isoname,evtype,new_file))
    return new_file

def getEnergycuts(ev):
    #open up event file and get header
    FT1=pyfits.open(ev)
    hdr=FT1[1].header
    #get number of keys in header and cycle through until we find the right one
    nkeys=hdr['NDSKEYS']
    for i in range(1,nkeys+1):
            if hdr['DSTYP%i'%i]=='ENERGY':
                    emin,emax=hdr['DSVAL%i'%i].split(':')
                    break
    #close the file and return
    FT1.close()
    return [float(emin),float(emax)]

def getROIcuts(ev):
    #open the file and get it's header
    FT1=pyfits.open(ev)
    hdr=FT1[1].header
    #get the number of keys in the header, cycle through until we find the right one
    nkeys=hdr['NDSKEYS']
    for i in range(1,nkeys+1):
        if hdr['DSTYP%i'%i]=='POS(RA,DEC)':
            #I'm not sure if this try statement is still necessary, doesn't hurt to leave it in though
            if "CIRCLE" in hdr['DSVAL%i'%i]:
                RA,DEC,RAD=hdr['DSVAL%i'%i].strip('CIRCLE()').split(',')
            else:
                RA,DEC,RAD=hdr['DSVAL%i'%i].strip('circle()').split(',')
            break
    #close the file and return
    FT1.close()
    return [float(RA),float(DEC),float(RAD)]


def type2bit(evt):
    if evt == "FRONT":
        return 1
    elif evt == "BACK":
        return 2
    elif "PSF" in evt:
        return 2**(int(evt[-1])+2)
    elif "EDISP" in evt:
        return 2**(int(evt[-1])+6)

def type2list(evtype):
    evtype_l = []
    if isinstance(evtype, basestring):
       if evtype not in ["FB","EDISP","PSF"]:
           if evtype in ALL_TYPES:
               evtype_l = [evtype]
           else:
               print "error"
       elif evtype=="FB":
           evtype_l = ["FRONT", "BACK"]
       else :
           evtype_l = ["%s%d"%(evtype,i) for i in range(4)]
    else:
        evtype_l = evtype
    return evtype_l


def wait_watching_stdout(ar, rc, dt=1, truncate=1000):
    while not ar.ready():
        rc.spin()
        stdouts = [ rc.metadata[msg_id]['stdout'] for msg_id in ar.msg_ids ]
        if not any(stdouts):
            continue
        # clear_output doesn't do much in terminal environments
        clear_output()
        print '-' * 30
        print "%.3fs elapsed" % ar.elapsed
        print ""
        for eid, stdout in zip(ar._targets, stdouts):
            if stdout:
                print "[ stdout %2i ]\n%s" % (eid, stdout[-truncate:])
        sys.stdout.flush()
        time.sleep(dt)
    
#@dview.remote(block=True)
def remote_gtselect():
    gtsel=GtApp('gtselect')
    outs=[]
    for type_name,type_bit  in evt_info:
        out = outbase + '_%s.fits'%type_name
        outs+=[out]
        
        if not os.access(out,os.F_OK):
            print 'Selecting %s %s events from %s, saving as %s.'%(evclass,type_name,ev,out)
            gtsel.run(infile=ev, outfile=out, ra=0.0,dec=0.0,\
                          rad=180., tmin=0., tmax=0., emin=0.,emax=0.,\
                          zmin=0.,zmax=0.,evclass=evclass,evtype=type_bit)
        else:
            print 'File %s already exists, skipping gtselect for %s.'%(out, ev)
    return

#@dview.remote(block=True)
def remote_gtbin():

    nebins = int(ceil((log10(ENE[1])-log10(ENE[0]))*10.))
    ebinalg = 'LOG'

    for type_name, type_bit  in evt_info:
        cc = os.path.join(workdir,'cc_%s_%s.fits'%(irfclass,type_name))
        #GTBIN:
        #INPUT : ENE, ROI, evfile, scfile
        #SCATTER OUTPUT: ccfiles
        #OTHER SCATTER? : binsz, to degrade maps with bad PSF or for BACK. How about ENE for EDISP?
        binsz=0.1
        nxpix=140
        nypix=140
        gtbin=GtApp('gtbin')
        if not os.access(cc,os.F_OK):
            gtbin.run(evfile=ev,scfile=scfile, outfile=cc, algorithm='CCUBE', \
                          ebinalg=ebinalg, emin=ENE[0], emax=ENE[1], enumbins=nebins, \
                          nxpix=nxpix, nypix=nypix, binsz=binsz, coordsys='CEL', \
                          xref=ROI[0], yref=ROI[1], axisrot=0, proj='AIT')
        else:
            print 'Counts Cube file %s already exists, skipping gtbin.'%cc

        #GTEXPCUBE2:
        #INPUT : ENE, ROI, ltc, ccfiles
        #SCATTER OUTPUT: ccfiles -> bdems
        #OTHER SCATTER? : binsz, to degrade maps with bad PSF or for BACK
        expcube_out = cc.replace('cc_','expcube_')
        #expcube_bin = int((ROI[2]+10.)*2./binsz)
        if not os.access(expcube_out,os.F_OK):
            bdem=GtApp('gtexpcube2')
            bdem.run(infile=ltc,cmap="none", outfile=expcube_out,\
                         irfs="P8R2_%s_V6"%irfclass, evtype=type_bit,\
                         nxpix=360, nypix=180,\
                         binsz=1, coordsys='CEL', xref=ROI[1],\
                         yref=ROI[2], axisrot=0, proj='AIT', emin=ENE[0],\
                         emax=ENE[1], enumbins=nebins, ebinalg=ebinalg)
        else:
            print '%s already exists, skipping gtexpcube2.'%expcube_out


def remote_srcmaps():
    #     #GTSRCMAPS:
    #     #INPUT : ltc, Expcube, ccfiles, diffuse models to scatter
    #     #SCATTER OUTPUT: ccfiles -> srcfiles
    #     #OTHER SCATTER? : binsz, to degrade maps with bad PSF or for BACK

    for type_name, type_bit  in evt_info:
        cc = os.path.join(workdir,'cc_%s_%s.fits'%(irfclass,type_name))
        expcube_out = cc.replace('cc_','expcube_')
        srcmap_out = cc.replace('cc_','srcmap_')

        if not os.access(srcmap_out,os.F_OK):
            srcmdl = build_xml(workdir, xml_basefile, type_name)
            gtsm=GtApp('gtsrcmaps')
            gtsm.run(expcube=ltc, cmap=cc, srcmdl=srcmdl, bexpmap=expcube_out, \
                             outfile=srcmap_out, irfs="CALDB")
        else:
            print 'Sourcemaps file %s already exists, skipping gtsrcmaps.'%srcmap_out


if __name__ == "__main__":
    import sys, yaml

    config_file = sys.argv[1]

    cfg = yaml.load(open(config_file))

    print "retrieving engines..."
    rc=Client()
    dview=rc[:]

    with dview.sync_imports():
        from GtApp import GtApp
        import os, sys, time, subprocess
        from math import ceil, log10

    print "launching the sequencer"

    irfclass = cfg['irfclass']
    workdir = cfg['workdir']
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    evfile = cfg['evfile']
    ltc = cfg['ltfile']
    evt_l = type2list(cfg['event_type'])
    evb_l = [type2bit(evt) for evt in evt_l]

    dview.push({
            'ltc':ltc,
            'irfclass':irfclass,
            'evclass':P8R2_CLASS_DICT[irfclass],
            'ev':evfile,
            'outbase':os.path.join(workdir,os.path.basename(evfile).split(".")[0]+'_'+irfclass),
            'workdir':workdir
            })
    dview.scatter('evt_info',zip(evt_l,evb_l))

    print "split events..."
    ar = dview.apply_async(remote_gtselect)
    wait_watching_stdout(ar,rc)
    split_ft1s = ar.get()
    if not ar.successful():
        print "ERROR in SELECT step"
    for am in ar.metadata:
        print am['stdout']

    #get ROI and Energy cuts from ev header
    ROI=getROIcuts(evfile)
    ENE=getEnergycuts(evfile)    
    dview.push({'ROI':ROI,'ENE':ENE,'scfile':cfg['scfile']})

    ar = dview.apply_async(remote_gtbin)
    wait_watching_stdout(ar, rc)
    ar.get()
    if not ar.successful():
        print "ERROR in BIN step"

    print "running srcmaps"
    dview.push({'xml_basefile':cfg['model'],
                'build_xml':build_xml
                })
    ar = dview.apply_async(remote_srcmaps)
    wait_watching_stdout(ar, rc)
    ar.get()
    if not ar.successful():
        print "ERROR in SRCMAP step"
