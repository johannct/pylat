import numpy as np
from constants import deg2rad

def angsep(x1_, y1_, x2_, y2_):
    """
    spherical angle separation, aka haversine formula
    input in degrees
    output is in degrees
    """
    deg2rad = pi/180.
    dlat = (y2_ - y1_)*deg2rad
    dlon = (x2_ - x1_)*deg2rad
    sina = np.sin(dlat/2.) * np.sin(dlat/2.) + \
        np.cos(y1_ * deg2rad) * np.cos(y2_ * deg2rad) * \
        np.sin(dlon/2.) * np.sin(dlon/2.)
    asina = 2 * np.arctan2(np.sqrt(sina), np.sqrt(1. - sina))
    return asina/deg2rad


def nearby_sources(analysis, RofI=360., RAcenter=None, DECcenter=None):
    """
    return the list of sources meeting the condition to be closer than RofI from position (RAcenter, Deccenter). If not provided, this center defaults to the center of the count map.
    """
    output = []
    for src in analysis.sourceNames():
        srcfuncs = analysis.model[src].funcs
        separation = -1.
        if RofI!=360. and 'Position' in srcfuncs.keys():
            position = srcfuncs['Position']
            ra = position.getParamValue("RA")
            dec = position.getParamValue("DEC")
            if RAcenter!=None and DECcenter!=None:
                separation = angsep(ra, dec, RAcenter, DECcenter)
            else:
                if hasattr(analysis,'components'):
                    c = analysis.components[0].logLike.countsMap()
                else:
                    c = analysis.logLike.countsMap()
                RAcenter  = c.crval1()
                DECcenter = c.crval2()
        if separation>0 and separation < RofI:
            output.append(src)
    return output
        
def set_free_params(analysis, srcName, params, free_map,verbose=False):
    """
    set the free/fixed flags of the params list to the corresponding flags in free_map, for source srcName. Length of free_map and params must match. 
    """
    if np.isscalar(params):
        params=[params]
        free_map=[free_map]
    spect = analysis.model[srcName].funcs['Spectrum']
    spect_type = spect.genericName()
    if len(params)!=len(free_map):
        print 'params and free_map lists must have the same lunch. Returning without action.'
    for i,pname in enumerate(params):
        try :
            spect.params[pname].setFree(free_map[i])
        except KeyError as e:
            print 'parameter %s in not part of %s parameters : '%(e,spect_type),spect.paramNames

    if verbose:
        print analysis[srcName]
    analysis.syncSrcParams()
    return
            

def modify_sources_roi_norm(analysis, srcName, free, RofI=360., RAcenter=None, DECcenter=None, calcNpred=False):
    """
    set free/fixed the norm parameter of source srcName.
    The source must be located closer than 
    RofI from position (RAcenter, DECcenter), if not None, or else 
    to the center of the countmaps. Default value for RofI is 360., which means that the condition is always satisfied.
    """
    srcfuncs = analysis.model[srcName].funcs
    separation = 0.
    if RofI!=360. and 'Position' in srcfuncs.keys():
        position = srcfuncs['Position']
        ra = position.getParamValue("RA")
        dec = position.getParamValue("DEC")
        if RAcenter!=None and DECcenter!=None:
            separation = angsep(ra, dec, RAcenter, DECcenter)
        else:
            c=analysis.logLike.countsMap()
            RAcenter  = c.crval1()
            DECcenter = c.crval2()
    if separation < RofI:
        analysis.normPar(srcName).setFree(free)
        #spectName = analysis[srcName].funcs['Spectrum'].genericName()

    analysis.logLike.syncSrcParams()        
    return

def get_model_map(analysis):
    """
    """
    mMap = pyLike.FloatVector()
    analysis.logLike.computeModelMap(mMap)
    a = np.asarray(mMap)
    cMap = analysis.logLike.countsMap()
    a.resize(len(cMap.energies())-1, cMap.naxis1(), cMap.naxis2())
    return a

