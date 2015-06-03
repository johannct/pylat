"""
This code subselect a list of sources from a Fermi LAT
catalog XML or DS9 region type file
based on a circular cut around a given position in the sky
"""
from xml.dom.minidom import parse
import numpy as np
from math import pi
import sys, os
from stlike_util import angsep

def rescale_normalisation(source, phase):
    """
    Look for the normalisation parameter and scale it by the
    phase factor
    """
    spectral_model = source.getElementsByTagName("spectrum")[0]
    norm_param = spectral_model.getElementsByTagName("parameter")[0]
    if norm_param.getAttribute("name") in ['Prefactor', 'Integral', 'norm']:
        val = float(norm_param.getAttribute("value"))
        newval = val*phase
        norm_param.setAttribute("value", str(newval))
    else:
        print "This should not happen:", \
            norm_param.getAttribute("name")
    return val, newval

def find_source_location(source):
    """
    look into the source element tag, which is a
    PointSource or DiffuseSource type, and extract
    the source ra and dec values
    """
    name = source.getAttribute("name")
    if source.getAttribute("type") == "PointSource":
        spatial_model = source.getElementsByTagName("spatialModel")[0]
        loc_params = spatial_model.getElementsByTagName("parameter")
        val1 = float(loc_params[0].getAttribute("value"))
        val2 = float(loc_params[1].getAttribute("value"))
    elif source.getAttribute("type") == "DiffuseSource":
        try:#the FGL entries have the RA/DEC info in the source element tag
            val1 = float(source.getAttribute("RA"))
            val2 = float(source.getAttribute("DEC"))
        except ValueError:
            return name, None, None
    return name, val1, val2

def skim_xml_file(filename, outfile, ra0, dec0, rad, phase):
    """
    This module deals with selection of sources from
    a ScienceTools gtlike model file
    """
    dom2 = parse(filename)
    main = dom2.childNodes[0] #element source_library
    for source in main.getElementsByTagName("source"):
        name, val1, val2 = find_source_location(source)
        if angsep(val1, val2, ra0, dec0) >= rad:
            main.removeChild(source)
        else:
            if phase != 1:
                val, newval = rescale_normalisation(source, phase)
                print name, val1, val2, angsep(val1, val2, ra0, dec0), \
                    "%g->%g" % (val, newval)
            else:
                print name, val1, val2, angsep(val1, val2, ra0, dec0)
        sys.stdout.flush()

    tempfile = "test.xml"
    fout = open(tempfile, "w")
    fout.write("<?xml version=\"1.0\" ?>\n")
    fout.write(dom2.saveXML(main))
    fout.close()
    #kludge to remove empty lines
    fout = open(tempfile)
    lines = fout.readlines()
    fout.close()
    while '\n' in lines:
        lines.remove('\n')
    fout = open(outfile, "w")
    fout.writelines(lines)
    fout.close()
    os.system('rm %s'%tempfile)

def skim_reg_file(fname, outfile, ra0, dec0, rad):
    """
    This module deals with selection of sources from
    a DS9 region file
    """
    fin_ = open(fname)
    out_ = open(outfile, "w")
    lines = fin_.readlines()
    for line in lines[1:]:
        if line[0] == "#":
            continue
        line2 = line.strip('fk5;point(   ')
        if line2[0]=="e":
            line2 = line2.strip('ellipse( ')
        ra_ = float(line2.split(",")[0])
        dec = float(line2.split(",")[1].split()[0].strip(')#'))
        id_ = line2.split()[-1].strip('text={}')
        if angsep(ra0, dec0, ra_, dec) <= rad:
            print ra_, dec, id_
            out_.write(line)
    fin_.close()
    out_.close()


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "usage: python %s X Y ROI infile outfile phase \
(infile is a reg or xml file, phase is the optional phase ratio to apply \
to all normalisations)" % sys.argv[0]
        exit()

    X_ = float(sys.argv[1])
    Y_ = float(sys.argv[2])
    ROI_ = float(sys.argv[3])
    IN_ = sys.argv[4]
    OUT_ = sys.argv[5]
    PHASE_ = 1
    if len(sys.argv) > 6:
        PHASE_ = float(sys.argv[6])

    try:
        skim_xml_file(IN_, OUT_, X_, Y_, ROI_, PHASE_)
    except Exception:
        skim_reg_file(IN_, OUT_, X_, Y_, ROI_)

